#pragma once
#include "animaPyramidalSymmetryConstrainedRegistrationBridge.h"

#include <itkTransformFactoryBase.h>
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkMultiResolutionPyramidImageFilter.h>
#include <itkImageRegistrationMethod.h>
#include <itkLinearInterpolateImageFunction.h>

#include <animaReadWriteFunctions.h>
#include <animaNewuoaOptimizer.h>
#include <animaMatrixOperations.h>
#include <animaResampleImageFilter.h>

#include <itkMeanSquaresImageToImageMetric.h>
#include <itkMutualInformationHistogramImageToImageMetric.h>
#include <itkNormalizedMutualInformationHistogramImageToImageMetric.h>
#include <itkCenteredTransformInitializer.h>

// ------------------------------

namespace anima
{

template <typename ScalarType>
PyramidalSymmetryConstrainedRegistrationBridge<ScalarType>::PyramidalSymmetryConstrainedRegistrationBridge()
{
    m_ReferenceImage = NULL;
    m_FloatingImage = NULL;

    m_OutputTransform = TransformType::New();
    m_OutputTransform->SetIdentity();

    m_outputTransformFile = "";

    m_InitialTransform = BaseTransformType::New();
    m_InitialTransform->SetIdentity();

    m_OutputImage = NULL;

    m_Metric = MutualInformation;
    m_OptimizerMaximumIterations = 100;

    m_SearchRadius = 2;
    m_SearchAngleRadius = 5;
    m_FinalRadius = 0.001;
    m_HistogramSize = 128;

    m_NumberOfPyramidLevels = 3;
    m_FastRegistration = false;
    this->SetNumberOfThreads(itk::MultiThreader::GetGlobalDefaultNumberOfThreads());
}

template <typename ScalarType>
PyramidalSymmetryConstrainedRegistrationBridge<ScalarType>::~PyramidalSymmetryConstrainedRegistrationBridge()
{
}

template <typename ScalarType>
void PyramidalSymmetryConstrainedRegistrationBridge<ScalarType>::Update()
{
    typedef typename itk::ImageRegistrationMethod<InputImageType, InputImageType> RegistrationType;
    this->SetupPyramids();

    typename TransformType::ParametersType initialParams(TransformType::ParametersDimension);
    for (unsigned int i = 0;i < TransformType::ParametersDimension;++i)
        initialParams[i] = 0;

    unsigned int indexAbsRefMax = 0;
    typename InputImageType::DirectionType dirRefMatrix = m_ReferenceImage->GetDirection();
    double valRefMax = std::abs(dirRefMatrix(0,0));
    for (unsigned int i = 1;i < InputImageType::ImageDimension;++i)
    {
        if (std::abs(dirRefMatrix(0,i)) > valRefMax)
        {
            valRefMax = std::abs(dirRefMatrix(0,i));
            indexAbsRefMax = i;
        }
    }
    
    OffsetType directionRefReal;
    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
        directionRefReal[i] = dirRefMatrix(i,indexAbsRefMax);

    m_OutputTransform->SetReferencePlaneNormal(directionRefReal);
    m_OutputTransform->SetParameters(initialParams);
    
    InputImageType::PointType centralPoint;
    itk::ContinuousIndex <ScalarType,InputImageType::ImageDimension> centralVoxIndex;
    
    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
        centralVoxIndex[i] = m_ReferenceImage->GetLargestPossibleRegion().GetSize()[i] / 2.0;
    
    m_ReferenceImage->TransformContinuousIndexToPhysicalPoint(centralVoxIndex,centralPoint);

    // Iterate over pyramid levels
    for (unsigned int i = 0;i < this->GetNumberOfPyramidLevels();++i)
    {
        std::cout << "Processing pyramid level " << i << std::endl;
        std::cout << "Image size: " << m_ReferencePyramid->GetOutput(i)->GetLargestPossibleRegion().GetSize() << std::endl;

        // Init matcher
        typename RegistrationType::Pointer reg = RegistrationType::New();

        reg->SetNumberOfThreads(this->GetNumberOfThreads());

        typedef anima::NewuoaOptimizer OptimizerType;

        typename OptimizerType::Pointer optimizer = OptimizerType::New();

        if (m_Metric == MeanSquares)
            optimizer->SetMaximize(false);
        else
            optimizer->SetMaximize(true);

        optimizer->SetRhoBegin(m_SearchRadius);
        optimizer->SetRhoEnd(m_FinalRadius);

        optimizer->SetNumberSamplingPoints(m_OutputTransform->GetNumberOfParameters() + 2);
        optimizer->SetMaximumIteration(m_OptimizerMaximumIterations);

        OptimizerType::ScalesType tmpScales (TransformType::ParametersDimension);
        tmpScales[0] = m_SearchRadius * 180.0 / (m_SearchAngleRadius * M_PI);
        
        unsigned int pos = 1;
        for (unsigned int j = 0;j < InputImageType::ImageDimension;++j)
        {
            if (j == indexAbsRefMax)
                continue;
            
            tmpScales[pos] = 1.0 / m_ReferencePyramid->GetOutput(i)->GetSpacing()[j];
            ++pos;
        }

        optimizer->SetScales(tmpScales);

        reg->SetOptimizer(optimizer);
        reg->SetTransform(m_OutputTransform);

        typedef itk::LinearInterpolateImageFunction <InputImageType,double> InterpolatorType;
        InterpolatorType::Pointer interpolator  = InterpolatorType::New();

        reg->SetInterpolator(interpolator);

        switch (m_Metric)
        {
            case MutualInformation:
            {
                typedef itk::MutualInformationHistogramImageToImageMetric < InputImageType,InputImageType > MetricType;
                typename MetricType::Pointer tmpMetric = MetricType::New();

                MetricType::HistogramType::SizeType histogramSize;
                histogramSize.SetSize(2);

                histogramSize[0] = m_HistogramSize;
                histogramSize[1] = m_HistogramSize;
                tmpMetric->SetHistogramSize( histogramSize );

                reg->SetMetric(tmpMetric);
                break;
            }

            case NormalizedMutualInformation:
            {
                typedef itk::NormalizedMutualInformationHistogramImageToImageMetric < InputImageType,InputImageType > MetricType;
                typename MetricType::Pointer tmpMetric = MetricType::New();

                MetricType::HistogramType::SizeType histogramSize;
                histogramSize.SetSize(2);

                histogramSize[0] = m_HistogramSize;
                histogramSize[1] = m_HistogramSize;
                tmpMetric->SetHistogramSize( histogramSize );

                reg->SetMetric(tmpMetric);
                break;
            }

            case MeanSquares:
            default:
            {
                typedef itk::MeanSquaresImageToImageMetric < InputImageType,InputImageType > MetricType;
                typename MetricType::Pointer tmpMetric = MetricType::New();
                reg->SetMetric(tmpMetric);
                break;
            }
        }

        reg->SetFixedImage(m_ReferencePyramid->GetOutput(i));
        reg->SetMovingImage(m_FloatingPyramid->GetOutput(i));

        if (m_FastRegistration)
        {
            // We can work in 2D since the rest should be perfectly ok
            InputImageRegionType workRegion = m_ReferencePyramid->GetOutput(i)->GetLargestPossibleRegion();
            InputImageRegionType::IndexType centralIndex;
            m_ReferencePyramid->GetOutput(i)->TransformPhysicalPointToIndex(centralPoint,centralIndex);

            unsigned int baseIndex = centralIndex[indexAbsRefMax];
            workRegion.SetIndex(indexAbsRefMax,baseIndex);
            workRegion.SetSize(indexAbsRefMax,1);

            reg->SetFixedImageRegion(workRegion);
        }
        else
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

        m_OutputTransform->SetParameters(reg->GetLastTransformParameters());
    }

    // Now compute the final transform
    typedef itk::Matrix <ScalarType,InputImageType::ImageDimension+1,InputImageType::ImageDimension+1> TransformMatrixType;
    TransformMatrixType refSymPlaneMatrix, initialMatrix, outputMatrix;
    refSymPlaneMatrix.SetIdentity();
    initialMatrix.SetIdentity();
    outputMatrix.SetIdentity();

    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
    {
        for (unsigned int j = 0;j < InputImageType::ImageDimension;++j)
        {
            refSymPlaneMatrix(i,j) = m_RefSymmetryTransform->GetMatrix()(i,j);
            outputMatrix(i,j) = m_OutputTransform->GetMatrix()(i,j);
            initialMatrix(i,j) = m_InitialTransform->GetMatrix()(i,j);
        }

        refSymPlaneMatrix(i,3) = m_RefSymmetryTransform->GetOffset()[i];
        outputMatrix(i,3) = m_OutputTransform->GetOffset()[i];
        initialMatrix(i,3) = m_InitialTransform->GetOffset()[i];
    }

    refSymPlaneMatrix = refSymPlaneMatrix.GetInverse();
    MatrixType tmpOutMatrix;
    OffsetType tmpOffset;

    outputMatrix = initialMatrix * outputMatrix * refSymPlaneMatrix;

    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
    {
        for (unsigned int j = 0;j < InputImageType::ImageDimension;++j)
            tmpOutMatrix(i,j) = outputMatrix(i,j);

        tmpOffset[i] = outputMatrix(i,3);
    }

    if (m_OutputRealignTransform.IsNull())
        m_OutputRealignTransform = BaseTransformType::New();

    m_OutputRealignTransform->SetMatrix(tmpOutMatrix);
    m_OutputRealignTransform->SetOffset(tmpOffset);

    typedef typename anima::ResampleImageFilter<InputImageType, InputImageType> ResampleFilterType;
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
}

template <typename ScalarType>
void PyramidalSymmetryConstrainedRegistrationBridge<ScalarType>::WriteOutputs()
{
    std::cout << "Writing output image to: " << m_resultFile << std::endl;

    anima::writeImage <InputImageType> (m_resultFile,m_OutputImage);

    if (m_outputTransformFile != "")
    {
        std::cout << "Writing output transform to: " << m_outputTransformFile << std::endl;
        itk::TransformFileWriter::Pointer writer = itk::TransformFileWriter::New();
        writer->SetInput(m_OutputRealignTransform);
        writer->SetFileName(m_outputTransformFile);
        writer->Update();
    }
}

template <typename ScalarType>
void PyramidalSymmetryConstrainedRegistrationBridge<ScalarType>::SetupPyramids()
{
    m_InitialTransform->SetIdentity();

    typedef typename itk::CenteredTransformInitializer<BaseTransformType, InputImageType, InputImageType> TransformInitializerType;
    typedef typename anima::ResampleImageFilter<InputImageType, InputImageType> ResampleFilterType;

    InputImagePointer initialReferenceImage;
    InputImagePointer initialFloatingImage;

    typename ResampleFilterType::Pointer tmpResample = ResampleFilterType::New();
    tmpResample->SetTransform(m_RefSymmetryTransform);
    tmpResample->SetInput(m_ReferenceImage);

    tmpResample->SetSize(m_ReferenceImage->GetLargestPossibleRegion().GetSize());
    tmpResample->SetOutputOrigin(m_ReferenceImage->GetOrigin());
    tmpResample->SetOutputSpacing(m_ReferenceImage->GetSpacing());
    tmpResample->SetOutputDirection(m_ReferenceImage->GetDirection());
    tmpResample->SetDefaultPixelValue(0);
    tmpResample->Update();

    initialReferenceImage = tmpResample->GetOutput();
    initialReferenceImage->DisconnectPipeline();

    // Tricky part: align the central planes of the two input images
    OffsetType directionRefReal, directionFloReal;

    //First use real direction to find the real X direction
    unsigned int indexAbsRefMax = 0;
    unsigned int indexAbsFloMax = 0;
    typename InputImageType::DirectionType dirRefMatrix = m_ReferenceImage->GetDirection();
    typename InputImageType::DirectionType dirFloMatrix = m_FloatingImage->GetDirection();
    double valRefMax = std::abs(dirRefMatrix(0,0));
    double valFloMax = std::abs(dirFloMatrix(0,0));
    for (unsigned int i = 1;i < InputImageType::ImageDimension;++i)
    {
        if (std::abs(dirRefMatrix(0,i)) > valRefMax)
        {
            valRefMax = std::abs(dirRefMatrix(0,i));
            indexAbsRefMax = i;
        }
        
        if (std::abs(dirFloMatrix(0,i)) > valFloMax)
        {
            valFloMax = std::abs(dirFloMatrix(0,i));
            indexAbsFloMax = i;
        }
    }

    // Now redo it with the real X-direction
    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
    {
        directionRefReal[i] = dirRefMatrix(i,indexAbsRefMax);
        directionFloReal[i] = dirFloMatrix(i,indexAbsFloMax);
    }

    typedef itk::Matrix <ScalarType,InputImageType::ImageDimension+1,InputImageType::ImageDimension+1> TransformMatrixType;
    MatrixType tmpMatrix = anima::GetRotationMatrixFromVectors(directionRefReal,directionFloReal);
    TransformMatrixType floRefMatrix;
    floRefMatrix.SetIdentity();
    
    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
    {
        for (unsigned int j = 0;j < InputImageType::ImageDimension;++j)
            floRefMatrix(i,j) = tmpMatrix(i,j);
    }

    itk::ContinuousIndex <ScalarType,InputImageType::ImageDimension> refImageCenter, floImageCenter;
    
    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
        refImageCenter[i] = m_ReferenceImage->GetLargestPossibleRegion().GetSize()[i] / 2.0;
    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
        floImageCenter[i] = m_FloatingImage->GetLargestPossibleRegion().GetSize()[i] / 2.0;

    typename InputImageType::PointType refCenter, floCenter;
    m_ReferenceImage->TransformContinuousIndexToPhysicalPoint(refImageCenter,refCenter);
    m_FloatingImage->TransformContinuousIndexToPhysicalPoint(floImageCenter,floCenter);

    TransformMatrixType refTranslationMatrix, floTranslationMatrix;
    refTranslationMatrix.SetIdentity();
    floTranslationMatrix.SetIdentity();
    
    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
    {
        refTranslationMatrix(i,3) = - refCenter[i];
        floTranslationMatrix(i,3) = floCenter[i];
    }

    floRefMatrix = floTranslationMatrix * floRefMatrix * refTranslationMatrix;
    
    TransformMatrixType floatingSymmetry;
    floatingSymmetry.SetIdentity();
    
    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
    {
        for (unsigned int j = 0;j < InputImageType::ImageDimension;++j)
            floatingSymmetry(i,j) = m_FloSymmetryTransform->GetMatrix()(i,j);

        floatingSymmetry(i,3) = m_FloSymmetryTransform->GetOffset()[i];
    }

    floRefMatrix = floatingSymmetry * floRefMatrix;
    m_OutputTransform->SetCenter(refCenter);

    MatrixType floatingMatrix;
    OffsetType floatingOffset;

    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
    {
        for (unsigned int j = 0;j < InputImageType::ImageDimension;++j)
            floatingMatrix(i,j) = floRefMatrix(i,j);

        floatingOffset[i] = floRefMatrix(i,3);
    }

    m_InitialTransform->SetMatrix(floatingMatrix);
    m_InitialTransform->SetOffset(floatingOffset);

    // Now, create pyramid
    m_ReferencePyramid = PyramidType::New();

    m_ReferencePyramid->SetInput(initialReferenceImage);
    m_ReferencePyramid->SetNumberOfLevels(m_NumberOfPyramidLevels);
    m_ReferencePyramid->Update();
    
    tmpResample = ResampleFilterType::New();
    tmpResample->SetTransform(m_InitialTransform);
    tmpResample->SetInput(m_FloatingImage);
    
    tmpResample->SetSize(m_ReferenceImage->GetLargestPossibleRegion().GetSize());
    tmpResample->SetOutputOrigin(m_ReferenceImage->GetOrigin());
    tmpResample->SetOutputSpacing(m_ReferenceImage->GetSpacing());
    tmpResample->SetOutputDirection(m_ReferenceImage->GetDirection());
    tmpResample->SetDefaultPixelValue(0);
    tmpResample->Update();
    
    initialFloatingImage = tmpResample->GetOutput();
    initialFloatingImage->DisconnectPipeline();

    // Create pyramid for floating image
    m_FloatingPyramid = PyramidType::New();

    m_FloatingPyramid->SetInput(initialFloatingImage);
    m_FloatingPyramid->SetNumberOfLevels(m_NumberOfPyramidLevels);
    m_FloatingPyramid->Update();
}

} // end of namespace
