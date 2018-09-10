#pragma once
#include "animaPyramidalSymmetryBridge.h"

#include <animaReadWriteFunctions.h>
#include <itkTransformFileWriter.h>
#include <itkMultiResolutionPyramidImageFilter.h>
#include <itkImageRegistrationMethod.h>
#include <itkLinearInterpolateImageFunction.h>
#include <animaNLOPTOptimizers.h>

#include <itkMeanSquaresImageToImageMetric.h>
#include <itkMutualInformationHistogramImageToImageMetric.h>
#include <itkImageMomentsCalculator.h>
#include <itkProgressReporter.h>

#include <animaVectorOperations.h>
#include <animaResampleImageFilter.h>

namespace anima
{

template <class PixelType, typename ScalarType>
void PyramidalSymmetryBridge<PixelType,ScalarType>::Update()
{
    typedef typename itk::ImageRegistrationMethod<OutputImageType, OutputImageType> RegistrationType;

    //progress management
    itk::ProgressReporter progress(this, 0, GetNumberOfPyramidLevels());

    if(m_progressCallback)
    {
        this->AddObserver ( itk::ProgressEvent(), m_progressCallback );
    }

    this->SetupPyramids();

    typename InputImageType::PointType centralPoint;

    typedef typename itk::ImageMomentsCalculator <InputImageType> ImageMomentsType;

    typename ImageMomentsType::Pointer momentsCalculator = ImageMomentsType::New();

    momentsCalculator->SetImage( m_ReferenceImage );
    momentsCalculator->Compute();

    itk::Vector <double,InputImageType::ImageDimension> centralVector = momentsCalculator->GetCenterOfGravity();

    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
        centralPoint[i] = centralVector[i];

    typename TransformType::ParametersType initialParams(TransformType::ParametersDimension);

    // Here should come test on orientation matrix -> find the true left/right axis
    for (unsigned int i = 0;i < TransformType::ParametersDimension;++i)
        initialParams[i] = 0;

    m_OutputTransform->SetParameters(initialParams);
    m_OutputTransform->SetRotationCenter(centralPoint);

    unsigned int dimension = m_OutputTransform->GetNumberOfParameters();
    itk::Array<double> lowerBounds(dimension);
    itk::Array<double> upperBounds(dimension);

    for (unsigned int i = 0;i < 2;++i)
    {
        lowerBounds[i] = - GetUpperBoundAngle();
        upperBounds[i] = GetUpperBoundAngle();
    }

    // Iterate over pyramid levels
    for (int i = 0;i < GetNumberOfPyramidLevels();++i)
    {
        std::cout << "Processing pyramid level " << i << std::endl;
        std::cout << "Image size: " << m_ReferencePyramid->GetOutput(i)->GetLargestPossibleRegion().GetSize() << std::endl;

        // Init matcher
        typename RegistrationType::Pointer reg = RegistrationType::New();

        reg->SetNumberOfWorkUnits(GetNumberOfWorkUnits());

        typedef anima::NLOPTOptimizers OptimizerType;
        typename OptimizerType::Pointer optimizer = OptimizerType::New();

        optimizer->SetAlgorithm(NLOPT_LN_BOBYQA);
        optimizer->SetXTolRel(1.0e-4);
        optimizer->SetFTolRel(1.0e-6);
        optimizer->SetMaxEval(GetOptimizerMaxIterations());
        optimizer->SetVectorStorageSize(2000);
        optimizer->SetMaximize(GetMetric() != MeanSquares);

        double meanSpacing = 0;
        for (unsigned int j = 0;j < InputImageType::ImageDimension;++j)
            meanSpacing += m_ReferencePyramid->GetOutput(i)->GetSpacing()[j];

        lowerBounds[2] = m_OutputTransform->GetParameters()[2] - meanSpacing * GetUpperBoundDistance();
        upperBounds[2] = m_OutputTransform->GetParameters()[2] + meanSpacing * GetUpperBoundDistance();

        optimizer->SetLowerBoundParameters(lowerBounds);
        optimizer->SetUpperBoundParameters(upperBounds);

        reg->SetOptimizer(optimizer);
        reg->SetTransform(m_OutputTransform);

        typedef itk::LinearInterpolateImageFunction <OutputImageType, double> InterpolatorType;
        typename InterpolatorType::Pointer interpolator = InterpolatorType::New();

        reg->SetInterpolator(interpolator);

        switch (GetMetric())
        {
            case MutualInformation:
            {
                typedef itk::MutualInformationHistogramImageToImageMetric < OutputImageType,OutputImageType > MetricType;
                typename MetricType::Pointer tmpMetric = MetricType::New();

                typename MetricType::HistogramType::SizeType histogramSize;
                histogramSize.SetSize(2);

                histogramSize[0] = GetHistogramSize();
                histogramSize[1] = GetHistogramSize();
                tmpMetric->SetHistogramSize( histogramSize );

                reg->SetMetric(tmpMetric);
                break;
            }
            case MeanSquares:
            default:
            {
                typedef itk::MeanSquaresImageToImageMetric < OutputImageType,OutputImageType > MetricType;
                typename MetricType::Pointer tmpMetric = MetricType::New();
                reg->SetMetric(tmpMetric);
                break;
            }
        }

        reg->SetFixedImage(m_ReferencePyramid->GetOutput(i));
        reg->SetMovingImage(m_FloatingPyramid->GetOutput(i));

        reg->SetFixedImageRegion(m_ReferencePyramid->GetOutput(i)->GetLargestPossibleRegion());
        reg->SetInitialTransformParameters(m_OutputTransform->GetParameters());

        try
        {
            reg->Update();
        }
        catch( itk::ExceptionObject & err )
        {
            std::cout << "ExceptionObject caught ! " << err << std::endl;
            throw err;
        }

        progress.CompletedPixel();
        m_OutputTransform->SetParameters(reg->GetLastTransformParameters());
    }

    // Now compute the transform to bring the image back onto its symmetry plane

    itk::ContinuousIndex <ScalarType,InputImageType::ImageDimension> imageCenter;

    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
        imageCenter[i] = (m_ReferenceImage->GetLargestPossibleRegion().GetSize()[i] - 1.0) / 2.0;

    typename InputImageType::PointType centerReal;
    m_ReferenceImage->TransformContinuousIndexToPhysicalPoint(imageCenter,centerReal);

    std::vector <double> directionVox(InputImageType::ImageDimension, 0);
    std::vector <double> directionReal(InputImageType::ImageDimension, 0);
    std::vector <double> directionSpherical(InputImageType::ImageDimension, 0);

    //First use real direction to find the real X direction
    unsigned int indexAbsMax = 0;
    typename InputImageType::DirectionType dirMatrix = m_ReferenceImage->GetDirection();
    double valMax = std::abs(dirMatrix(0,0));
    for (unsigned int i = 1;i < InputImageType::ImageDimension;++i)
    {
        if (std::abs(dirMatrix(0,i)) > valMax)
        {
            valMax = std::abs(dirMatrix(0,i));
            indexAbsMax = i;
        }
    }

    // Now redo it with the real X-direction
    directionVox[indexAbsMax] = 1;

    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
        directionReal[i] = dirMatrix(i,indexAbsMax);

    anima::TransformCartesianToSphericalCoordinates(directionReal,directionSpherical);

    initialParams.Fill(0);
    initialParams[0] = M_PI / 2.0 - directionSpherical[0];
    initialParams[1] = directionSpherical[1];
    initialParams[2] = 0;

    this->ComputeRealignTransform(centralVector,centerReal,initialParams);

    // Compute output image
    typedef typename anima::ResampleImageFilter<InputImageType, OutputImageType> ResampleFilterType;
    typename ResampleFilterType::Pointer tmpResample = ResampleFilterType::New();
    tmpResample->SetTransform(m_OutputRealignTransform);
    tmpResample->SetInput(m_FloatingImage);

    tmpResample->SetSize(m_ReferenceImage->GetLargestPossibleRegion().GetSize());
    tmpResample->SetOutputOrigin(m_ReferenceImage->GetOrigin());
    tmpResample->SetOutputSpacing(m_ReferenceImage->GetSpacing());
    tmpResample->SetOutputDirection(m_ReferenceImage->GetDirection());
    tmpResample->SetDefaultPixelValue(0);
    tmpResample->Update();

    m_OutputImage = tmpResample->GetOutput();
    m_OutputImage->DisconnectPipeline();
}

template <class PixelType, typename ScalarType>
void PyramidalSymmetryBridge<PixelType,ScalarType>::ComputeRealignTransform(typename itk::Vector <double,InputImageType::ImageDimension> centralPoint,
                                                                            typename InputImageType::PointType &centerReal, ParametersType &imageParams)
{
    TransformPointer imageMidPlaneTrsf = TransformType::New();
    imageMidPlaneTrsf->SetParameters(imageParams);
    imageMidPlaneTrsf->SetRotationCenter(centerReal);

    typedef itk::Matrix <ScalarType,InputImageType::ImageDimension+1,InputImageType::ImageDimension+1> HomogeneousMatrixType;

    HomogeneousMatrixType midPlaneMatrix;
    midPlaneMatrix.SetIdentity();

    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
    {
        midPlaneMatrix(i,InputImageType::ImageDimension) = imageMidPlaneTrsf->GetOffset()[i];

        for (unsigned int j = 0;j < InputImageType::ImageDimension;++j)
            midPlaneMatrix(i,j) = imageMidPlaneTrsf->GetMatrix()(i,j);
    }

    HomogeneousMatrixType outputTrsfMatrix;
    outputTrsfMatrix.SetIdentity();

    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
    {
        outputTrsfMatrix(i,InputImageType::ImageDimension) = m_OutputTransform->GetOffset()[i];

        for (unsigned int j = 0;j < InputImageType::ImageDimension;++j)
            outputTrsfMatrix(i,j) = m_OutputTransform->GetMatrix()(i,j);
    }

    HomogeneousMatrixType matrixComposition = midPlaneMatrix * outputTrsfMatrix;

    double theta_x, theta_y, theta_z;
    double sinus, cosinus, phi, f;

    // note : le vecteur rotation obtenu est toujours dans le plan
    // au milieu de l'image, donc une des 3 composantes est forcement nulle
    // ce qui donne l'idee d'une autre parametrisation du plan de symetrie

    theta_x = - matrixComposition[1][2] + matrixComposition[2][1]; //X = R(2,3) - R(3,2)
    theta_y = - matrixComposition[2][0] + matrixComposition[0][2]; //Y = R(3,1) - R(1,3)
    theta_z = - matrixComposition[0][1] + matrixComposition[1][0]; //Z = R(1,2) - R(2,1)

    sinus = std::sqrt(theta_x*theta_x + theta_y*theta_y + theta_z*theta_z);

    if (std::abs(sinus) > 1e-9)
    {
        cosinus = matrixComposition[0][0] + matrixComposition[1][1] + matrixComposition[2][2] - 1;
        phi = std::atan(sinus / cosinus);
        f = phi / sinus;
        theta_x = theta_x * f;
        theta_y = theta_y * f;
        theta_z = theta_z * f;
    }
    else
    {
        phi = 0;

        theta_x = 0;
        theta_y = 0;
        theta_z = 0;
    }

    // le vecteur rotation de r est celui de R divisé par 2

    phi = phi / 2.0;
    theta_x = theta_x / 2.0;
    theta_y = theta_y / 2.0;
    theta_z = theta_z / 2.0;

    // formules : http://www.iau-sofa.rl.ac.uk/2003_0429/sofa/rv2m.for

    // calcul du vecteur unitaire ayant la même direction que le vecteur rotation

    double u_x, u_y, u_z;

    sinus = std::sqrt(theta_x*theta_x + theta_y*theta_y + theta_z*theta_z);

    if (std::abs(sinus) > 1e-9)
    {
        u_x = theta_x / sinus;
        u_y = theta_y / sinus;
        u_z = theta_z / sinus;
    }

    else {
        u_x = 0;
        u_y = 0;
        u_z = 0;
    }

    // calcul de r

    double c, s, t;

    c = std::cos(phi);
    s = std::sin(phi);
    t = 1 - c;

    HomogeneousMatrixType sqrtMatrix;
    sqrtMatrix.SetIdentity();

    sqrtMatrix[0][0] = t*u_x*u_x + c;
    sqrtMatrix[0][1] = t*u_x*u_y - s*u_z;
    sqrtMatrix[0][2] = t*u_x*u_z + s*u_y;
    sqrtMatrix[1][0] = t*u_x*u_y + s*u_z;
    sqrtMatrix[1][1] = t*u_y*u_y + c;
    sqrtMatrix[1][2] = t*u_y*u_z - s*u_x;
    sqrtMatrix[2][0] = t*u_x*u_z - s*u_y;
    sqrtMatrix[2][1] = t*u_y*u_z + s*u_x;
    sqrtMatrix[2][2] = t*u_z*u_z + c;

    // calcul de transfo5 = r+I

    MatrixType tmpMatrix;

    for (unsigned int i=0; i<3; i++)
        for (unsigned int j=0; j<3; j++)
            tmpMatrix[i][j] = sqrtMatrix[i][j];

    for (unsigned int i=0; i<3; i++)
        tmpMatrix[i][i]++;

    // calcul de transfo6 = (r+I)^{-1}

    tmpMatrix = tmpMatrix.GetInverse();

    // calcul de t = (r+I)^{-1}.T
    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
        for (unsigned int j = 0;j < InputImageType::ImageDimension;++j)
            sqrtMatrix[i][3] += tmpMatrix[i][j] * matrixComposition[j][3];

    itk::Vector <double,InputImageType::ImageDimension+1> centralPointHom, centralPointTr;
    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
        centralPointHom[i] = centralPoint[i];
    centralPointHom[InputImageType::ImageDimension] = 1;

    centralPointTr = sqrtMatrix * centralPointHom;

    sqrtMatrix = sqrtMatrix.GetInverse();

    HomogeneousMatrixType trToCenterMatrix;
    trToCenterMatrix.SetIdentity();

    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
        trToCenterMatrix(i,3) = centralPointTr[i] - centerReal[i];

    sqrtMatrix = sqrtMatrix * trToCenterMatrix;

    MatrixType trsfMatrix;
    OffsetType offsetVector;
    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
    {
        offsetVector[i] = sqrtMatrix(i,3);
        for (unsigned int j=0; j<3; j++)
            trsfMatrix(i,j) = sqrtMatrix(i,j);
    }

    if (m_OutputRealignTransform.IsNull())
        m_OutputRealignTransform = BaseTransformType::New();

    m_OutputRealignTransform->SetMatrix(trsfMatrix);
    m_OutputRealignTransform->SetOffset(offsetVector);
}

template <class PixelType, typename ScalarType>
void PyramidalSymmetryBridge<PixelType,ScalarType>::WriteOutputs()
{
    SaveResultFile();

    SaveRealignTransformFile();

    SaveTransformFile();
}


template <class PixelType, typename ScalarType>
void PyramidalSymmetryBridge<PixelType,ScalarType>::SaveResultFile()
{
    if (GetResultfile() != "")
    {
        std::cout << "Writing output image to: " << GetResultfile() << std::endl;
        anima::writeImage <InputImageType> (GetResultfile(),m_OutputImage);
    }
}


template <class PixelType, typename ScalarType>
void PyramidalSymmetryBridge<PixelType,ScalarType>::SaveTransformFile()
{
    if (GetOutputTransformFile() != "")
    {
        std::cout << "Writing output transform to: " << GetOutputTransformFile() << std::endl;
        itk::TransformFileWriter::Pointer writer = itk::TransformFileWriter::New();

        // SymmetryPlaneTransforms should not be written as is, this loses information
        typename BaseTransformType::Pointer tmpTrsf = BaseTransformType::New();
        tmpTrsf->SetMatrix(m_OutputTransform->GetMatrix());
        tmpTrsf->SetOffset(m_OutputTransform->GetOffset());

        writer->SetInput(tmpTrsf);
        writer->SetFileName(GetOutputTransformFile());
        writer->Update();
    }
}

template <class PixelType, typename ScalarType>
void PyramidalSymmetryBridge<PixelType,ScalarType>::SaveRealignTransformFile()
{
    if (GetOutputRealignTransformFile() != "")
    {
        std::cout << "Writing output realign transform to: " << GetOutputRealignTransformFile() << std::endl;
        itk::TransformFileWriter::Pointer writer = itk::TransformFileWriter::New();
        writer->SetInput(m_OutputRealignTransform);
        writer->SetFileName(GetOutputRealignTransformFile());
        writer->Update();
    }
}

template <class PixelType, typename ScalarType>
void PyramidalSymmetryBridge<PixelType,ScalarType>::SetupPyramids()
{
    m_ReferencePyramid = PyramidType::New();

    m_ReferencePyramid->SetInput(m_ReferenceImage);
    m_ReferencePyramid->SetNumberOfLevels(GetNumberOfPyramidLevels());

    m_ReferencePyramid->Update();

    // Create pyramid for floating image
    m_FloatingPyramid = PyramidType::New();

    m_FloatingPyramid->SetInput(m_FloatingImage);
    m_FloatingPyramid->SetNumberOfLevels(GetNumberOfPyramidLevels());
    m_FloatingPyramid->Update();
}

} // end of namespace anima
