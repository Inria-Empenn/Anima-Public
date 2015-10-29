#pragma once
#include "animaPyramidalDistortionCorrectionBlockMatchingBridge.h"
#include <itkImageMomentsCalculator.h>

#include <animaDistortionCorrectionBMRegistrationMethod.h>
#include <animaDistortionCorrectionBlockMatcher.h>
#include <animaResampleImageFilter.h>

#include <itkVectorResampleImageFilter.h>

#include <animaVelocityUtils.h>
#include <animaReadWriteFunctions.h>

namespace anima
{

template <unsigned int ImageDimension>
PyramidalDistortionCorrectionBlockMatchingBridge<ImageDimension>::PyramidalDistortionCorrectionBlockMatchingBridge()
{
    m_BackwardImage = NULL;
    m_ForwardImage = NULL;
    
    m_InitialTransform = NULL;
    
    m_OutputTransform = DisplacementFieldTransformType::New();
    m_OutputTransform->SetIdentity();
    
    m_outputTransformFile = "";
    
    m_OutputImage = NULL;
    
    m_TransformDirection = 0;
    m_BlockSize = 5;
    m_BlockSpacing = 2;
    m_StDevThreshold = 5;
    m_MaximumIterations = 10;
    m_OptimizerMaximumIterations = 100;
    m_SearchRadius = 2;
    m_SearchScaleRadius = 0.1;
    m_SearchSkewRadius = 0.1;
    m_FinalRadius = 0.001;
    m_TranlateUpperBound = 50;
    m_ScaleUpperBound = std::log(5.0);
    m_SkewUpperBound = std::tan(M_PI / 3.0);
    m_Agregator = Baloo;
    m_TransformKind = DirectionScaleSkew;
    m_Metric = SquaredCorrelation;
    m_WeightedAgregation = false;
    m_ExtrapolationSigma = 3;
    m_ElasticSigma = 3;
    m_OutlierSigma = 3;
    m_MEstimateConvergenceThreshold = 0.01;
    m_NeighborhoodApproximation = 2.5;
    m_UseTransformationDam = true;
    m_DamDistance = 2.5;
    m_NumberOfPyramidLevels = 3;
    m_LastPyramidLevel = 0;
    m_PercentageKept = 0.8;
    this->SetNumberOfThreads(itk::MultiThreader::GetGlobalDefaultNumberOfThreads());
}

template <unsigned int ImageDimension>
void
PyramidalDistortionCorrectionBlockMatchingBridge<ImageDimension>::
SetInitialTransformField(VectorFieldType *field)
{
    if (!m_InitialTransform)
        m_InitialTransform = DisplacementFieldTransformType::New();

    m_InitialTransform->SetParametersAsVectorField(field);
}

template <unsigned int ImageDimension>
PyramidalDistortionCorrectionBlockMatchingBridge<ImageDimension>::~PyramidalDistortionCorrectionBlockMatchingBridge()
{
}

template <unsigned int ImageDimension>
void 
PyramidalDistortionCorrectionBlockMatchingBridge<ImageDimension>::Update()
{
    typedef typename anima::DistortionCorrectionBMRegistrationMethod<InputImageType> BlockMatchRegistrationType;

    this->SetupPyramids();

    // Iterate over pyramid levels
    for (unsigned int i = 0;i < m_BackwardPyramid->GetNumberOfLevels();++i)
    {
        if (i + m_LastPyramidLevel >= m_BackwardPyramid->GetNumberOfLevels())
            continue;

        typename InputImageType::Pointer backwardImage = m_BackwardPyramid->GetOutput(i);
        backwardImage->DisconnectPipeline();
        
        typename InputImageType::Pointer forwardImage = m_ForwardPyramid->GetOutput(i);
        forwardImage->DisconnectPipeline();
        
        // Update fields to match the current resolution
        if (m_OutputTransform->GetParametersAsVectorField() != NULL)
        {
            typedef itk::VectorResampleImageFilter<VectorFieldType,VectorFieldType> VectorResampleFilterType;
            typedef typename VectorResampleFilterType::Pointer VectorResampleFilterPointer;
            
            AffineTransformPointer tmpIdentity = AffineTransformType::New();
            tmpIdentity->SetIdentity();
            
            VectorResampleFilterPointer tmpResample = VectorResampleFilterType::New();
            tmpResample->SetTransform(tmpIdentity);
            tmpResample->SetInput(m_OutputTransform->GetParametersAsVectorField());
            
            tmpResample->SetSize(backwardImage->GetLargestPossibleRegion().GetSize());
            tmpResample->SetOutputOrigin(backwardImage->GetOrigin());
            tmpResample->SetOutputSpacing(backwardImage->GetSpacing());
            tmpResample->SetOutputDirection(backwardImage->GetDirection());
            
            tmpResample->Update();
            
            VectorFieldType *tmpOut = tmpResample->GetOutput();
            m_OutputTransform->SetParametersAsVectorField(tmpOut);
            tmpOut->DisconnectPipeline();
        }
        
        DisplacementFieldTransformPointer initialTransform;
        
        if (m_InitialTransform)
        {
            typedef itk::VectorResampleImageFilter<VectorFieldType,VectorFieldType> VectorResampleFilterType;
            typedef typename VectorResampleFilterType::Pointer VectorResampleFilterPointer;
            
            AffineTransformPointer tmpIdentity = AffineTransformType::New();
            tmpIdentity->SetIdentity();
            
            VectorResampleFilterPointer tmpResample = VectorResampleFilterType::New();
            tmpResample->SetTransform(tmpIdentity);
            tmpResample->SetInput(m_InitialTransform->GetParametersAsVectorField());
            
            tmpResample->SetSize(backwardImage->GetLargestPossibleRegion().GetSize());
            tmpResample->SetOutputOrigin(backwardImage->GetOrigin());
            tmpResample->SetOutputSpacing(backwardImage->GetSpacing());
            tmpResample->SetOutputDirection(backwardImage->GetDirection());
            
            tmpResample->Update();
            
            VectorFieldType *tmpOut = tmpResample->GetOutput();
            initialTransform = DisplacementFieldTransformType::New();
            initialTransform->SetParametersAsVectorField(tmpOut);
            tmpOut->DisconnectPipeline();
        }
        
        std::cout << "Processing pyramid level " << i << std::endl;
        std::cout << "Image size: " << backwardImage->GetLargestPossibleRegion().GetSize() << std::endl;
        
        double meanSpacing = 0;
        for (unsigned int j = 0;j < ImageDimension;++j)
            meanSpacing += backwardImage->GetSpacing()[j];
        meanSpacing /= ImageDimension;
        
        // Init matcher
        typename BlockMatchRegistrationType::Pointer bmreg = BlockMatchRegistrationType::New();

        bmreg->SetNumberOfThreads(this->GetNumberOfThreads());
        
        typedef anima::ResampleImageFilter<InputImageType, InputImageType,
                                         typename BaseAgregatorType::ScalarType> ResampleFilterType;

        typename ResampleFilterType::Pointer refResampler = ResampleFilterType::New();
        refResampler->SetSize(forwardImage->GetLargestPossibleRegion().GetSize());
        refResampler->SetOutputOrigin(forwardImage->GetOrigin());
        refResampler->SetOutputSpacing(forwardImage->GetSpacing());
        refResampler->SetOutputDirection(forwardImage->GetDirection());
        refResampler->SetDefaultPixelValue(0);
        refResampler->SetNumberOfThreads(GetNumberOfThreads());
        refResampler->SetScaleIntensitiesWithJacobian(true);
        bmreg->SetReferenceImageResampler(refResampler);

        typename ResampleFilterType::Pointer movingResampler = ResampleFilterType::New();
        movingResampler->SetSize(backwardImage->GetLargestPossibleRegion().GetSize());
        movingResampler->SetOutputOrigin(backwardImage->GetOrigin());
        movingResampler->SetOutputSpacing(backwardImage->GetSpacing());
        movingResampler->SetOutputDirection(backwardImage->GetDirection());
        movingResampler->SetDefaultPixelValue(0);
        movingResampler->SetNumberOfThreads(GetNumberOfThreads());
        movingResampler->SetScaleIntensitiesWithJacobian(true);
        bmreg->SetMovingImageResampler(movingResampler);

        // Init matcher
        typedef anima::DistortionCorrectionBlockMatcher <InputImageType> BlockMatcherType;

        BlockMatcherType *mainMatcher = new BlockMatcherType;
        mainMatcher->SetBlockPercentageKept(GetPercentageKept());
        mainMatcher->SetBlockSize(GetBlockSize());
        mainMatcher->SetBlockSpacing(GetBlockSpacing());
        mainMatcher->SetBlockVarianceThreshold(GetStDevThreshold() * GetStDevThreshold());
        mainMatcher->SetUseTransformationDam(m_UseTransformationDam);
        mainMatcher->SetDamDistance(m_DamDistance * m_ExtrapolationSigma);

        bmreg->SetBlockMatcher(mainMatcher);

        bmreg->SetFixedImage(backwardImage);
        bmreg->SetMovingImage(forwardImage);
        
        // Init agregator mean shift parameters
        BaseAgregatorType* agregPtr = NULL;
        
        if (m_Agregator == MSmoother)
        {
            MEstimateAgregatorType *agreg = new MEstimateAgregatorType;
            agreg->SetExtrapolationSigma(m_ExtrapolationSigma * meanSpacing);
            agreg->SetOutlierRejectionSigma(m_OutlierSigma);
            agreg->SetOutputTransformType(BaseAgregatorType::SVF);
            
            agreg->SetNumberOfThreads(this->GetNumberOfThreads());
            agreg->SetGeometryInformation(backwardImage.GetPointer());
            
            agreg->SetNeighborhoodHalfSize((unsigned int)floor(m_ExtrapolationSigma * m_NeighborhoodApproximation));
            agreg->SetDistanceBoundary(m_ExtrapolationSigma * meanSpacing * m_NeighborhoodApproximation);
            agreg->SetMEstimateConvergenceThreshold(m_MEstimateConvergenceThreshold);

            agregPtr = agreg;
        }
        else
        {
            BalooAgregatorType *agreg = new BalooAgregatorType;
            agreg->SetExtrapolationSigma(m_ExtrapolationSigma * meanSpacing);
            agreg->SetOutlierRejectionSigma(m_OutlierSigma);
            agreg->SetOutputTransformType(BaseAgregatorType::SVF);
            
            agreg->SetNumberOfThreads(this->GetNumberOfThreads());
            agreg->SetGeometryInformation(backwardImage.GetPointer());
            
            agregPtr = agreg;
        }

        bmreg->SetAgregator(agregPtr);

        mainMatcher->SetBlockTransformType((typename BlockMatcherType::TransformDefinition) m_TransformKind);
        mainMatcher->SetSimilarityType((typename BlockMatcherType::SimilarityDefinition) m_Metric);
        mainMatcher->SetOptimizerType(BlockMatcherType::Bobyqa);

        bmreg->SetSVFElasticRegSigma(m_ElasticSigma * meanSpacing);

        bmreg->SetMaximumIterations(m_MaximumIterations);
        mainMatcher->SetOptimizerMaximumIterations(m_OptimizerMaximumIterations);
        mainMatcher->SetTransformDirection(m_TransformDirection);
        
        if (initialTransform)
            bmreg->SetInitialTransform(initialTransform.GetPointer());
        
        bmreg->SetCurrentTransform(m_OutputTransform.GetPointer());

        mainMatcher->SetSearchRadius(m_SearchRadius);
        mainMatcher->SetSearchScaleRadius(m_SearchScaleRadius);
        mainMatcher->SetSearchSkewRadius(m_SearchSkewRadius);
        mainMatcher->SetFinalRadius(m_FinalRadius);
        mainMatcher->SetTranslateMax(m_TranlateUpperBound);
        mainMatcher->SetScaleMax(m_ScaleUpperBound);
        mainMatcher->SetSkewMax(m_SkewUpperBound);
        mainMatcher->SetNumberOfThreads(GetNumberOfThreads());

        bmreg->Update();

        const DisplacementFieldTransformType *resTrsf = dynamic_cast <const DisplacementFieldTransformType *> (bmreg->GetOutput()->Get());
        m_OutputTransform->SetParametersAsVectorField(resTrsf->GetParametersAsVectorField());
    }
    
    if (m_LastPyramidLevel != 0)
    {
        // Resample output transform to go back to full resolution
        typedef itk::VectorResampleImageFilter<VectorFieldType,VectorFieldType> VectorResampleFilterType;
        typedef typename VectorResampleFilterType::Pointer VectorResampleFilterPointer;
        
        AffineTransformPointer tmpIdentity = AffineTransformType::New();
        tmpIdentity->SetIdentity();
        
        VectorResampleFilterPointer tmpResample = VectorResampleFilterType::New();
        tmpResample->SetTransform(tmpIdentity);
        tmpResample->SetInput(m_OutputTransform->GetParametersAsVectorField());
        
        tmpResample->SetSize(m_BackwardImage->GetLargestPossibleRegion().GetSize());
        tmpResample->SetOutputOrigin(m_BackwardImage->GetOrigin());
        tmpResample->SetOutputSpacing(m_BackwardImage->GetSpacing());
        tmpResample->SetOutputDirection(m_BackwardImage->GetDirection());
        
        tmpResample->Update();
        
        VectorFieldType *tmpOut = tmpResample->GetOutput();
        m_OutputTransform->SetParametersAsVectorField(tmpOut);
        tmpOut->DisconnectPipeline();
    }
    
    typedef typename anima::ResampleImageFilter<InputImageType, InputImageType,
            typename BaseAgregatorType::ScalarType> ResampleFilterType;
    
    DisplacementFieldTransformPointer oppositeTransform = DisplacementFieldTransformType::New();
    
    typedef itk::MultiplyImageFilter <VectorFieldType,itk::Image <float, ImageDimension>, VectorFieldType> MultiplyFilterType;
    typename MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
    multiplyFilter->SetInput(m_OutputTransform->GetParametersAsVectorField());
    multiplyFilter->SetConstant(-1.0);
    multiplyFilter->SetNumberOfThreads(this->GetNumberOfThreads());
    
    multiplyFilter->Update();

    oppositeTransform->SetParametersAsVectorField(multiplyFilter->GetOutput());
    
    if (m_InitialTransform)
        anima::composeDistortionCorrections<typename BaseAgregatorType::ScalarType,ImageDimension>
                (m_InitialTransform,m_OutputTransform,oppositeTransform,this->GetNumberOfThreads());
    
    typename ResampleFilterType::Pointer tmpResampleFloating = ResampleFilterType::New();
    tmpResampleFloating->SetTransform(m_OutputTransform);
    tmpResampleFloating->SetInput(m_ForwardImage);
    
    tmpResampleFloating->SetSize(m_BackwardImage->GetLargestPossibleRegion().GetSize());
    tmpResampleFloating->SetOutputOrigin(m_BackwardImage->GetOrigin());
    tmpResampleFloating->SetOutputSpacing(m_BackwardImage->GetSpacing());
    tmpResampleFloating->SetOutputDirection(m_BackwardImage->GetDirection());
    tmpResampleFloating->SetDefaultPixelValue(0);
    tmpResampleFloating->SetScaleIntensitiesWithJacobian(true);
    tmpResampleFloating->SetNumberOfThreads(this->GetNumberOfThreads());
    tmpResampleFloating->Update();
    
    typename ResampleFilterType::Pointer tmpResampleReference = ResampleFilterType::New();
    tmpResampleReference->SetTransform(oppositeTransform);
    tmpResampleReference->SetInput(m_BackwardImage);
    
    tmpResampleReference->SetSize(m_BackwardImage->GetLargestPossibleRegion().GetSize());
    tmpResampleReference->SetOutputOrigin(m_BackwardImage->GetOrigin());
    tmpResampleReference->SetOutputSpacing(m_BackwardImage->GetSpacing());
    tmpResampleReference->SetOutputDirection(m_BackwardImage->GetDirection());
    tmpResampleReference->SetDefaultPixelValue(0);
    tmpResampleReference->SetScaleIntensitiesWithJacobian(true);
    tmpResampleReference->SetNumberOfThreads(this->GetNumberOfThreads());
    
    tmpResampleReference->Update();
    
    typedef itk::AddImageFilter <InputImageType,InputImageType,InputImageType> AddFilterType;
    typename AddFilterType::Pointer addFilter = AddFilterType::New();
    addFilter->SetInput1(tmpResampleFloating->GetOutput());
    addFilter->SetInput2(tmpResampleReference->GetOutput());
    addFilter->SetNumberOfThreads(this->GetNumberOfThreads());
    
    addFilter->Update();
    
    typedef itk::MultiplyImageFilter <InputImageType,itk::Image <float, ImageDimension>,InputImageType> MultiplyScalarFilterType;
    typename MultiplyScalarFilterType::Pointer multiplyScalarFilter = MultiplyScalarFilterType::New();
    multiplyScalarFilter->SetInput(addFilter->GetOutput());
    multiplyScalarFilter->SetConstant(0.5);
    multiplyScalarFilter->SetNumberOfThreads(this->GetNumberOfThreads());
    multiplyScalarFilter->InPlaceOn();
    
    multiplyScalarFilter->Update();
    
    m_OutputImage = multiplyScalarFilter->GetOutput();
    m_OutputImage->DisconnectPipeline();
}

template <unsigned int ImageDimension>
void 
PyramidalDistortionCorrectionBlockMatchingBridge<ImageDimension>::WriteOutputs()
{
    std::cout << "Writing output image to: " << m_resultFile << std::endl;
    anima::writeImage<InputImageType>(m_resultFile,m_OutputImage);

    if (m_outputTransformFile != "")
    {
        std::cout << "Writing output transform to: " << m_outputTransformFile << std::endl;
        anima::writeImage<VectorFieldType>(m_outputTransformFile, const_cast <VectorFieldType *> (m_OutputTransform->GetParametersAsVectorField()));
    }
}

template <unsigned int ImageDimension>
void 
PyramidalDistortionCorrectionBlockMatchingBridge<ImageDimension>::SetupPyramids()
{
    // Create pyramid here, check images actually are of the same size.
    m_BackwardPyramid = PyramidType::New();

    m_BackwardPyramid->SetInput(m_BackwardImage);
    m_BackwardPyramid->SetNumberOfLevels(m_NumberOfPyramidLevels);
    m_BackwardPyramid->SetNumberOfThreads(this->GetNumberOfThreads());

    typedef typename anima::ResampleImageFilter<InputImageType, InputImageType,
            typename BaseAgregatorType::ScalarType> ResampleFilterType;

    typename ResampleFilterType::Pointer backwardResampler = ResampleFilterType::New();
    m_BackwardPyramid->SetImageResampler(backwardResampler);

    m_BackwardPyramid->Update();

    // Create pyramid for floating image
    m_ForwardPyramid = PyramidType::New();

    m_ForwardPyramid->SetInput(m_ForwardImage);
    m_ForwardPyramid->SetNumberOfLevels(m_NumberOfPyramidLevels);
    m_ForwardPyramid->SetNumberOfThreads(this->GetNumberOfThreads());

    typename ResampleFilterType::Pointer forwardResampler = ResampleFilterType::New();
    m_ForwardPyramid->SetImageResampler(forwardResampler);

    m_ForwardPyramid->Update();
}

} // end namespace anima
