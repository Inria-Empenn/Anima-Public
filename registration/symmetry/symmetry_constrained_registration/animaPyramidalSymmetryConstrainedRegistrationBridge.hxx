#pragma once

#include "animaPyramidalSymmetryConstrainedRegistrationBridge.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTransformFactoryBase.h>
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkMultiResolutionPyramidImageFilter.h>
#include <itkImageRegistrationMethod.h>
#include <itkLinearInterpolateImageFunction.h>
#include <animaNewuoaOptimizer.h>
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

    m_OutputTransform->SetReferenceSymmetryPlanes(m_RefSymmetryTransform,m_FloSymmetryTransform);
    m_OutputTransform->SetParameters(initialParams);

    typedef typename itk::ImageMomentsCalculator <InputImageType> ImageMomentsType;
    typename ImageMomentsType::Pointer momentsCalculator = ImageMomentsType::New();

    momentsCalculator->SetImage(m_ReferenceImage);
    momentsCalculator->Compute();

    itk::Vector <double,InputImageType::ImageDimension> centralVector = momentsCalculator->GetCenterOfGravity();
    typename InputImageType::PointType centralPoint;
    typename InputImageType::IndexType centralIndex;
    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
        centralPoint[i] = centralVector[i];

    // Iterate over pyramid levels
    for (unsigned int i = 0;i < m_NumberOfPyramidLevels;++i)
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

        OptimizerType::ScalesType tmpScales( TransformType::ParametersDimension );
        tmpScales[0] = m_SearchRadius * 180.0 / (m_SearchAngleRadius * M_PI);
        tmpScales[1] = 1.0 / m_ReferencePyramid->GetOutput(i)->GetSpacing()[1];
        tmpScales[2] = 1.0 / m_ReferencePyramid->GetOutput(i)->GetSpacing()[2];

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
            m_ReferencePyramid->GetOutput(i)->TransformPhysicalPointToIndex(centralPoint,centralIndex);

            unsigned int baseIndex = centralIndex[0];
            workRegion.SetIndex(0,baseIndex);
            workRegion.SetSize(0,1);

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

        std::cout << "Converged after " << optimizer->GetCurrentIteration() << std::endl;

        m_OutputTransform->SetParameters(reg->GetLastTransformParameters());
    }

    // Now compute the final transform

    itk::Matrix <ScalarType,InputImageType::ImageDimension+1,InputImageType::ImageDimension+1> initialMatrix, outputMatrix;
    initialMatrix.SetIdentity();
    outputMatrix.SetIdentity();

    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
    {
        for (unsigned int j = 0;j < InputImageType::ImageDimension;++j)
        {
            outputMatrix(i,j) = m_OutputTransform->GetMatrix()(i,j);
            initialMatrix(i,j) = m_InitialTransform->GetMatrix()(i,j);
        }

        outputMatrix(i,3) = m_OutputTransform->GetOffset()[i];
        initialMatrix(i,3) = m_InitialTransform->GetOffset()[i];
    }

    MatrixType tmpOutMatrix;
    OffsetType tmpOffset;

    outputMatrix = initialMatrix * outputMatrix;

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

    typedef typename itk::ResampleImageFilter<InputImageType, InputImageType> ResampleFilterType;
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

    typename itk::ImageFileWriter <InputImageType>::Pointer imageWriter = itk::ImageFileWriter <InputImageType>::New();
    imageWriter->SetUseCompression(true);
    imageWriter->SetInput(m_OutputImage);
    imageWriter->SetFileName(m_resultFile);

    imageWriter->Update();

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
    typedef typename itk::ResampleImageFilter<InputImageType, InputImageType> ResampleFilterType;

    InputImagePointer initialFloatingImage = m_FloatingImage;

    typename TransformInitializerType::Pointer initializer = TransformInitializerType::New();
    initializer->SetTransform(m_InitialTransform);
    initializer->SetFixedImage(m_ReferenceImage);
    initializer->SetMovingImage(m_FloatingImage);
    initializer->MomentsOn();
    initializer->InitializeTransform();

    typename ResampleFilterType::Pointer tmpResample = ResampleFilterType::New();
    tmpResample->SetTransform(m_InitialTransform);
    tmpResample->SetInput(m_FloatingImage);

    tmpResample->SetSize(m_ReferenceImage->GetLargestPossibleRegion().GetSize());
    tmpResample->SetOutputOrigin(m_ReferenceImage->GetOrigin());
    tmpResample->SetOutputSpacing(m_ReferenceImage->GetSpacing());
    tmpResample->SetOutputDirection(m_ReferenceImage->GetDirection());
    tmpResample->SetDefaultPixelValue(0);
    tmpResample->Update();

    initialFloatingImage = tmpResample->GetOutput();

    // Now the ugly part, basically we need the images to have the same sizes for pyramids but then we need to update the floating
    // symmetry to make it work

    itk::Matrix <ScalarType,InputImageType::ImageDimension+1,InputImageType::ImageDimension+1> initialMatrix, tmpMatrix, floatingSymmetry;
    initialMatrix.SetIdentity();
    floatingSymmetry.SetIdentity();
    tmpMatrix.SetIdentity();

    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
    {
        for (unsigned int j = 0;j < InputImageType::ImageDimension;++j)
        {
            initialMatrix(i,j) = m_InitialTransform->GetMatrix()(i,j);
            floatingSymmetry(i,j) = m_FloSymmetryTransform->GetMatrix()(i,j);
        }

        initialMatrix(i,3) = m_InitialTransform->GetOffset()[i];
        floatingSymmetry(i,3) = m_FloSymmetryTransform->GetOffset()[i];
    }

    initialMatrix = initialMatrix.GetInverse();

    floatingSymmetry = initialMatrix * floatingSymmetry;

    itk::ContinuousIndex <ScalarType,InputImageType::ImageDimension> refImageCenter, floImageCenter;

    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
        refImageCenter[i] = m_ReferenceImage->GetLargestPossibleRegion().GetSize()[i] / 2.0;
    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
        floImageCenter[i] = m_FloatingImage->GetLargestPossibleRegion().GetSize()[i] / 2.0;

    typename InputImageType::PointType refCenter, floCenter;
    m_ReferenceImage->TransformContinuousIndexToPhysicalPoint(refImageCenter,refCenter);
    m_FloatingImage->TransformContinuousIndexToPhysicalPoint(floImageCenter,floCenter);

    m_OutputTransform->SetRotationCenter(refCenter);

    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
        tmpMatrix(i,3) = floCenter[i] - refCenter[i];

    floatingSymmetry = floatingSymmetry * tmpMatrix;

    MatrixType floatingMatrix;
    OffsetType floatingOffset;

    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
    {
        for (unsigned int j = 0;j < InputImageType::ImageDimension;++j)
            floatingMatrix(i,j) = floatingSymmetry(i,j);

        floatingOffset[i] = floatingSymmetry(i,3);
    }

    m_FloSymmetryTransform->SetMatrix(floatingMatrix);
    m_FloSymmetryTransform->SetOffset(floatingOffset);

    // Now, create pyramid
    typedef itk::ResampleImageFilter<InputImageType, InputImageType> ResampleFilterType;

    m_ReferencePyramid = PyramidType::New();

    m_ReferencePyramid->SetInput(m_ReferenceImage);
    m_ReferencePyramid->SetNumberOfLevels(m_NumberOfPyramidLevels);
    m_ReferencePyramid->Update();

    // Create pyramid for floating image
    m_FloatingPyramid = PyramidType::New();

    m_FloatingPyramid->SetInput(initialFloatingImage);
    m_FloatingPyramid->SetNumberOfLevels(m_NumberOfPyramidLevels);
    m_FloatingPyramid->Update();
}

} // end of namespace
