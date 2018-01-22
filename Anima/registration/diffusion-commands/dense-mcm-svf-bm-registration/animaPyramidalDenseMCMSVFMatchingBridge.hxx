#pragma once
#include "animaPyramidalDenseMCMSVFMatchingBridge.h"

#include <animaReadWriteFunctions.h>
#include <animaMCMFileWriter.h>

#include <animaAsymmetricBMRegistrationMethod.h>
#include <animaSymmetricBMRegistrationMethod.h>
#include <animaKissingSymmetricBMRegistrationMethod.h>

#include <itkVectorResampleImageFilter.h>

#include <animaVelocityUtils.h>
#include <animaMCMResampleImageFilter.h>

namespace anima
{

template <unsigned int ImageDimension>
PyramidalDenseMCMSVFMatchingBridge<ImageDimension>::PyramidalDenseMCMSVFMatchingBridge()
{
    m_ReferenceImage = NULL;
    m_FloatingImage = NULL;

    m_OutputTransform = BaseTransformType::New();
    m_OutputTransform->SetIdentity();

    m_outputTransformFile = "";

    m_OutputImage = NULL;

    m_BlockSize = 5;
    m_BlockSpacing = 2;
    m_StDevThreshold = 5;

    m_SymmetryType = Asymmetric;
    m_MetricOrientation = FiniteStrain;
    m_FiniteStrainImageReorientation = true;
    m_Transform = Translation;
    m_Metric = MCMOneToOneBasicMeanSquares;
    m_Optimizer = Bobyqa;

    m_MaximumIterations = 10;
    m_MinimalTransformError = 0.01;
    m_OptimizerMaximumIterations = 100;
    m_SearchRadius = 2;
    m_SearchAngleRadius = 5;
    m_SearchScaleRadius = 0.1;
    m_FinalRadius = 0.001;
    m_StepSize = 1;
    m_TranslateUpperBound = 50;
    m_AngleUpperBound = 180;
    m_ScaleUpperBound = 3;
    m_Agregator = Baloo;
    m_ExtrapolationSigma = 3;
    m_ElasticSigma = 3;
    m_OutlierSigma = 3;
    m_MEstimateConvergenceThreshold = 0.01;
    m_NeighborhoodApproximation = 2.5;
    m_BCHCompositionOrder = 1;
    m_ExponentiationOrder = 1;
    m_UseTransformationDam = true;
    m_DamDistance = 2.5;
    m_NumberOfPyramidLevels = 3;
    m_LastPyramidLevel = 0;
    m_PercentageKept = 0.8;
    this->SetNumberOfThreads(itk::MultiThreader::GetGlobalDefaultNumberOfThreads());
}

template <unsigned int ImageDimension>
PyramidalDenseMCMSVFMatchingBridge<ImageDimension>::~PyramidalDenseMCMSVFMatchingBridge()
{
}

template <unsigned int ImageDimension>
typename PyramidalDenseMCMSVFMatchingBridge<ImageDimension>::InterpolatorType *
PyramidalDenseMCMSVFMatchingBridge<ImageDimension>
::CreateInterpolator(InputImageType *image)
{
    typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
    interpolator->SetReferenceOutputModel(image->GetDescriptionModel());

    interpolator->Register();
    return interpolator;
}

template <unsigned int ImageDimension>
typename PyramidalDenseMCMSVFMatchingBridge<ImageDimension>::BlockMatcherType *
PyramidalDenseMCMSVFMatchingBridge<ImageDimension>
::CreateBlockMatcher()
{
    BlockMatcherType *matcher = new BlockMatcherType;
    return matcher;
}

template <unsigned int ImageDimension>
void
PyramidalDenseMCMSVFMatchingBridge<ImageDimension>::Update()
{
    this->SetupPyramids();

    // Iterate over pyramid levels
    for (unsigned int i = 0;i < m_ReferencePyramid->GetNumberOfLevels();++i)
    {
        if (i + m_LastPyramidLevel >= m_ReferencePyramid->GetNumberOfLevels())
            continue;

        typename InputImageType::Pointer refImage = m_ReferencePyramid->GetOutput(i);
        refImage->DisconnectPipeline();

        typename InputImageType::Pointer floImage = m_FloatingPyramid->GetOutput(i);
        floImage->DisconnectPipeline();

        // Update field to match the current resolution
        if (m_OutputTransform->GetParametersAsVectorField() != NULL)
        {
            typedef itk::VectorResampleImageFilter<VelocityFieldType,VelocityFieldType> VectorResampleFilterType;
            typedef typename VectorResampleFilterType::Pointer VectorResampleFilterPointer;

            AffineTransformPointer tmpIdentity = AffineTransformType::New();
            tmpIdentity->SetIdentity();

            VectorResampleFilterPointer tmpResample = VectorResampleFilterType::New();
            tmpResample->SetTransform(tmpIdentity);
            tmpResample->SetInput(m_OutputTransform->GetParametersAsVectorField());

            tmpResample->SetSize(refImage->GetLargestPossibleRegion().GetSize());
            tmpResample->SetOutputOrigin(refImage->GetOrigin());
            tmpResample->SetOutputSpacing(refImage->GetSpacing());
            tmpResample->SetOutputDirection(refImage->GetDirection());

            tmpResample->Update();

            VelocityFieldType *tmpOut = tmpResample->GetOutput();
            m_OutputTransform->SetParametersAsVectorField(tmpOut);
            tmpOut->DisconnectPipeline();
        }

        std::cout << "Processing pyramid level " << i << std::endl;
        std::cout << "Image size: " << refImage->GetLargestPossibleRegion().GetSize() << std::endl;

        double meanSpacing = 0;
        for (unsigned int j = 0;j < ImageDimension;++j)
            meanSpacing += refImage->GetSpacing()[j];
        meanSpacing /= ImageDimension;

        // Init agregator mean shift parameters
        BaseAgregatorType* agregPtr = NULL;

        if (m_Agregator == MSmoother)
        {
            MEstimateAgregatorType *agreg = new MEstimateAgregatorType;
            agreg->SetExtrapolationSigma(m_ExtrapolationSigma * meanSpacing);
            agreg->SetOutlierRejectionSigma(m_OutlierSigma);
            agreg->SetOutputTransformType(BaseAgregatorType::SVF);

            if (this->GetNumberOfThreads() != 0)
                agreg->SetNumberOfThreads(this->GetNumberOfThreads());

            agreg->SetGeometryInformation(refImage.GetPointer());

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

            if (this->GetNumberOfThreads() != 0)
                agreg->SetNumberOfThreads(this->GetNumberOfThreads());

            agreg->SetGeometryInformation(refImage.GetPointer());

            agregPtr = agreg;
        }

        // Init matcher
        BlockMatcherType *mainMatcher = this->CreateBlockMatcher();

        BlockMatcherType *reverseMatcher = 0;
        mainMatcher->SetBlockPercentageKept(GetPercentageKept());
        mainMatcher->SetBlockSize(GetBlockSize());
        mainMatcher->SetBlockSpacing(GetBlockSpacing());
        mainMatcher->SetBlockVarianceThreshold(GetStDevThreshold() * GetStDevThreshold());
        mainMatcher->SetUseTransformationDam(m_UseTransformationDam);
        mainMatcher->SetDamDistance(m_DamDistance * meanSpacing / 2.0);
        mainMatcher->SetGradientDirections(m_GradientDirections);
        mainMatcher->SetBValues(m_BValues);

        switch (m_SymmetryType)
        {
            case Asymmetric:
            {
                typedef typename anima::AsymmetricBMRegistrationMethod <InputImageType> BlockMatchRegistrationType;
                m_bmreg = BlockMatchRegistrationType::New();
                break;
            }

            case Symmetric:
            {
                typedef typename anima::SymmetricBMRegistrationMethod <InputImageType> BlockMatchRegistrationType;
                typename BlockMatchRegistrationType::Pointer tmpReg = BlockMatchRegistrationType::New();
                reverseMatcher = this->CreateBlockMatcher();
                reverseMatcher->SetBlockPercentageKept(GetPercentageKept());
                reverseMatcher->SetBlockSize(GetBlockSize());
                reverseMatcher->SetBlockSpacing(GetBlockSpacing());
                reverseMatcher->SetBlockVarianceThreshold(GetStDevThreshold() * GetStDevThreshold());
                reverseMatcher->SetUseTransformationDam(m_UseTransformationDam);
                reverseMatcher->SetDamDistance(m_DamDistance * meanSpacing / 2.0);
                reverseMatcher->SetGradientDirections(m_GradientDirections);
                reverseMatcher->SetBValues(m_BValues);

                tmpReg->SetReverseBlockMatcher(reverseMatcher);
                m_bmreg = tmpReg;
                break;
            }

            case Kissing:
            {
                typedef typename anima::KissingSymmetricBMRegistrationMethod <InputImageType> BlockMatchRegistrationType;
                m_bmreg = BlockMatchRegistrationType::New();
                break;
            }
        }

        m_bmreg->SetBlockMatcher(mainMatcher);
        m_bmreg->SetAgregator(agregPtr);
        m_bmreg->SetBCHCompositionOrder(m_BCHCompositionOrder);
        m_bmreg->SetExponentiationOrder(m_ExponentiationOrder);

        if (this->GetNumberOfThreads() != 0)
            m_bmreg->SetNumberOfThreads(this->GetNumberOfThreads());

        m_bmreg->SetFixedImage(refImage);
        m_bmreg->SetMovingImage(floImage);

        m_bmreg->SetSVFElasticRegSigma(m_ElasticSigma * meanSpacing);

        typedef typename anima::MCMResampleImageFilter<InputImageType,typename BaseAgregatorType::ScalarType> ResampleFilterType;

        typename InterpolatorType::Pointer interpolator = this->CreateInterpolator(floImage);

        typename ResampleFilterType::Pointer refResampler = ResampleFilterType::New();
        refResampler->SetOutputLargestPossibleRegion(floImage->GetLargestPossibleRegion());
        refResampler->SetOutputOrigin(floImage->GetOrigin());
        refResampler->SetOutputSpacing(floImage->GetSpacing());
        refResampler->SetOutputDirection(floImage->GetDirection());
        refResampler->SetNumberOfThreads(GetNumberOfThreads());
        refResampler->SetReferenceOutputModel(floImage->GetDescriptionModel());
        refResampler->SetFiniteStrainReorientation(this->GetFiniteStrainImageReorientation());
        refResampler->SetInterpolator(interpolator.GetPointer());
        m_bmreg->SetReferenceImageResampler(refResampler);

        interpolator = this->CreateInterpolator(refImage);

        typename ResampleFilterType::Pointer movingResampler = ResampleFilterType::New();
        movingResampler->SetOutputLargestPossibleRegion(refImage->GetLargestPossibleRegion());
        movingResampler->SetOutputOrigin(refImage->GetOrigin());
        movingResampler->SetOutputSpacing(refImage->GetSpacing());
        movingResampler->SetOutputDirection(refImage->GetDirection());
        movingResampler->SetNumberOfThreads(GetNumberOfThreads());
        movingResampler->SetReferenceOutputModel(refImage->GetDescriptionModel());
        movingResampler->SetFiniteStrainReorientation(this->GetFiniteStrainImageReorientation());
        movingResampler->SetInterpolator(interpolator.GetPointer());
        m_bmreg->SetMovingImageResampler(movingResampler);

        switch (GetTransform())
        {
            case Translation:
                mainMatcher->SetBlockTransformType(BlockMatcherType::Superclass::Translation);
                if (reverseMatcher)
                    reverseMatcher->SetBlockTransformType(BlockMatcherType::Superclass::Translation);
                break;
            case Rigid:
                mainMatcher->SetBlockTransformType(BlockMatcherType::Superclass::Rigid);
                if (reverseMatcher)
                    reverseMatcher->SetBlockTransformType(BlockMatcherType::Superclass::Rigid);
                break;
            case Affine:
            default:
                mainMatcher->SetBlockTransformType(BlockMatcherType::Superclass::Affine);
                if (reverseMatcher)
                    reverseMatcher->SetBlockTransformType(BlockMatcherType::Superclass::Affine);
                break;
        }

        switch (GetOptimizer())
        {
            case Exhaustive:
                mainMatcher->SetOptimizerType(BlockMatcherType::Exhaustive);
                if (reverseMatcher)
                    reverseMatcher->SetOptimizerType(BlockMatcherType::Exhaustive);
                break;

            case Bobyqa:
            default:
                mainMatcher->SetOptimizerType(BlockMatcherType::Bobyqa);
                if (reverseMatcher)
                    reverseMatcher->SetOptimizerType(BlockMatcherType::Bobyqa);
                break;
        }

        switch (m_Metric)
        {
            case MCMBasicMeanSquares:
                mainMatcher->SetSimilarityType(BlockMatcherType::MCMBasicMeanSquares);
                if (reverseMatcher)
                    reverseMatcher->SetSimilarityType(BlockMatcherType::MCMBasicMeanSquares);
                break;
            case MCMOneToOneBasicMeanSquares:
                mainMatcher->SetSimilarityType(BlockMatcherType::MCMOneToOneBasicMeanSquares);
                if (reverseMatcher)
                    reverseMatcher->SetSimilarityType(BlockMatcherType::MCMOneToOneBasicMeanSquares);
                break;
            case MCMCorrelation:
                mainMatcher->SetSimilarityType(BlockMatcherType::MCMCorrelation);
                if (reverseMatcher)
                    reverseMatcher->SetSimilarityType(BlockMatcherType::MCMCorrelation);
                break;
            case MTCorrelation:
                mainMatcher->SetSimilarityType(BlockMatcherType::MTCorrelation);
                if (reverseMatcher)
                    reverseMatcher->SetSimilarityType(BlockMatcherType::MTCorrelation);
                break;
            case MCMMeanSquares:
            default:
                mainMatcher->SetSimilarityType(BlockMatcherType::MCMMeanSquares);
                if (reverseMatcher)
                    reverseMatcher->SetSimilarityType(BlockMatcherType::MCMMeanSquares);
                break;
        }

        switch (m_MetricOrientation)
        {
            case None:
                mainMatcher->SetModelRotationType(BlockMatcherType::None);
                if (reverseMatcher)
                    reverseMatcher->SetModelRotationType(BlockMatcherType::None);
                break;
            case PPD:
                mainMatcher->SetModelRotationType(BlockMatcherType::PPD);
                if (reverseMatcher)
                    reverseMatcher->SetModelRotationType(BlockMatcherType::PPD);
                break;
            case FiniteStrain:
            default:
                mainMatcher->SetModelRotationType(BlockMatcherType::FiniteStrain);
                if (reverseMatcher)
                    reverseMatcher->SetModelRotationType(BlockMatcherType::FiniteStrain);
                break;
        }

        m_bmreg->SetMaximumIterations(m_MaximumIterations);
        m_bmreg->SetMinimalTransformError(m_MinimalTransformError);
        m_bmreg->SetInitialTransform(m_OutputTransform.GetPointer());

        mainMatcher->SetOptimizerMaximumIterations(m_OptimizerMaximumIterations);

        double sr = m_SearchRadius;
        mainMatcher->SetSearchRadius(sr);

        double sar = m_SearchAngleRadius;
        mainMatcher->SetSearchAngleRadius(sar);

        double scr = m_SearchScaleRadius;
        mainMatcher->SetSearchScaleRadius(scr);

        double fr = m_FinalRadius;
        mainMatcher->SetFinalRadius(fr);

        double ss = m_StepSize;
        mainMatcher->SetStepSize(ss);

        double tub = m_TranslateUpperBound;
        mainMatcher->SetTranslateMax(tub);

        double aub = m_AngleUpperBound;
        mainMatcher->SetAngleMax(aub);

        double scub = m_ScaleUpperBound;
        mainMatcher->SetScaleMax(scub);

        if (reverseMatcher)
        {
            reverseMatcher->SetNumberOfThreads(GetNumberOfThreads());
            reverseMatcher->SetOptimizerMaximumIterations(GetOptimizerMaximumIterations());

            reverseMatcher->SetSearchRadius(sr);
            reverseMatcher->SetSearchAngleRadius(sar);
            reverseMatcher->SetSearchScaleRadius(scr);
            reverseMatcher->SetFinalRadius(fr);
            reverseMatcher->SetStepSize(ss);
            reverseMatcher->SetTranslateMax(tub);
            reverseMatcher->SetAngleMax(aub);
            reverseMatcher->SetScaleMax(scub);
        }

        try
        {
            m_bmreg->Update();
        }
        catch (itk::ExceptionObject & err)
        {
            std::cout << "Exception: " << err << std::endl;
            exit(-1);
        }

        const BaseTransformType *resTrsf = dynamic_cast <const BaseTransformType *> (m_bmreg->GetOutput()->Get());
        m_OutputTransform->SetParametersAsVectorField(resTrsf->GetParametersAsVectorField());

        delete mainMatcher;
        if (reverseMatcher)
            delete reverseMatcher;
    }
    
    if (m_SymmetryType == Kissing)
    {
        VelocityFieldType *finalTrsfField = const_cast <VelocityFieldType *> (m_OutputTransform->GetParametersAsVectorField());
        typedef itk::MultiplyImageFilter <VelocityFieldType,itk::Image <float,ImageDimension>, VelocityFieldType> MultiplyFilterType;

        typename MultiplyFilterType::Pointer fieldMultiplier = MultiplyFilterType::New();
        fieldMultiplier->SetInput(finalTrsfField);
        fieldMultiplier->SetConstant(2.0);
        fieldMultiplier->SetNumberOfThreads(GetNumberOfThreads());
        fieldMultiplier->InPlaceOn();

        fieldMultiplier->Update();

        VelocityFieldType *outputField = fieldMultiplier->GetOutput();
        m_OutputTransform->SetParametersAsVectorField(fieldMultiplier->GetOutput());
        outputField->DisconnectPipeline();
    }

    DisplacementFieldTransformPointer outputDispTrsf = DisplacementFieldTransformType::New();
    anima::GetSVFExponential(m_OutputTransform.GetPointer(), outputDispTrsf.GetPointer(), m_ExponentiationOrder, GetNumberOfThreads(), false);

    typedef typename anima::MCMResampleImageFilter<InputImageType,typename BaseAgregatorType::ScalarType> ResampleFilterType;
    typename ResampleFilterType::Pointer tmpResample = ResampleFilterType::New();

    typename InterpolatorType::Pointer interpolator = this->CreateInterpolator(m_ReferenceImage);

    typedef itk::Transform<typename BaseAgregatorType::ScalarType,ImageDimension,ImageDimension> BaseTransformType;
    typename BaseTransformType::Pointer baseTrsf = outputDispTrsf.GetPointer();
    tmpResample->SetTransform(baseTrsf);
    tmpResample->SetFiniteStrainReorientation(this->GetFiniteStrainImageReorientation());
    tmpResample->SetInput(m_FloatingImage);

    if (this->GetNumberOfThreads() != 0)
        tmpResample->SetNumberOfThreads(this->GetNumberOfThreads());

    tmpResample->SetOutputLargestPossibleRegion(m_ReferenceImage->GetLargestPossibleRegion());
    tmpResample->SetOutputOrigin(m_ReferenceImage->GetOrigin());
    tmpResample->SetOutputSpacing(m_ReferenceImage->GetSpacing());
    tmpResample->SetOutputDirection(m_ReferenceImage->GetDirection());
    tmpResample->SetReferenceOutputModel(m_ReferenceImage->GetDescriptionModel());
    tmpResample->SetInterpolator(interpolator.GetPointer());
    tmpResample->Update();

    m_OutputImage = tmpResample->GetOutput();
    m_OutputImage->DisconnectPipeline();
}

template <unsigned int ImageDimension>
void
PyramidalDenseMCMSVFMatchingBridge<ImageDimension>::WriteOutputs()
{
    std::cout << "Writing output image to: " << m_resultFile << std::endl;
    anima::MCMFileWriter <float,ImageDimension> mcmWriter;
    mcmWriter.SetInputImage(m_OutputImage);
    mcmWriter.SetFileName(m_resultFile);

    mcmWriter.Update();

    if (m_outputTransformFile != "")
    {
        std::cout << "Writing output SVF to: " << m_outputTransformFile << std::endl;
        anima::writeImage <VelocityFieldType> (m_outputTransformFile,
                                               const_cast <VelocityFieldType *> (m_OutputTransform->GetParametersAsVectorField()));
    }
}

template <unsigned int ImageDimension>
void
PyramidalDenseMCMSVFMatchingBridge<ImageDimension>::SetupPyramids()
{
    // Create pyramid here, check images actually are of the same size.
    m_ReferencePyramid = PyramidType::New();

    m_ReferencePyramid->SetInput(m_ReferenceImage);
    m_ReferencePyramid->SetNumberOfLevels(m_NumberOfPyramidLevels);

    if (this->GetNumberOfThreads() != 0)
        m_ReferencePyramid->SetNumberOfThreads(this->GetNumberOfThreads());

    typedef typename anima::MCMResampleImageFilter<InputImageType,typename BaseAgregatorType::ScalarType> ResampleFilterType;
    typename InterpolatorType::Pointer interpolator = this->CreateInterpolator(m_ReferenceImage);

    typename ResampleFilterType::Pointer refResampler = ResampleFilterType::New();
    refResampler->SetReferenceOutputModel(m_ReferenceImage->GetDescriptionModel());
    refResampler->SetFiniteStrainReorientation(this->GetFiniteStrainImageReorientation());
    refResampler->SetInterpolator(interpolator.GetPointer());
    m_ReferencePyramid->SetImageResampler(refResampler);

    m_ReferencePyramid->Update();

    // Create pyramid for floating image
    m_FloatingPyramid = PyramidType::New();

    m_FloatingPyramid->SetInput(m_FloatingImage);
    m_FloatingPyramid->SetNumberOfLevels(m_NumberOfPyramidLevels);

    if (this->GetNumberOfThreads() != 0)
        m_FloatingPyramid->SetNumberOfThreads(this->GetNumberOfThreads());

    typename ResampleFilterType::Pointer floResampler = ResampleFilterType::New();
    interpolator = this->CreateInterpolator(m_FloatingImage);

    floResampler->SetReferenceOutputModel(m_FloatingImage->GetDescriptionModel());
    floResampler->SetFiniteStrainReorientation(this->GetFiniteStrainImageReorientation());
    floResampler->SetInterpolator(interpolator.GetPointer());
    m_FloatingPyramid->SetImageResampler(floResampler);

    m_FloatingPyramid->Update();
}

} // end of namespace anima
