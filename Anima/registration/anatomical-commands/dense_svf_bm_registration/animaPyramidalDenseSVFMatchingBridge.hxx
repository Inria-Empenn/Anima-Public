#pragma once
#include "animaPyramidalDenseSVFMatchingBridge.h"

#include <itkResampleImageFilter.h>
#include <itkMinimumMaximumImageFilter.h>

#include <animaReadWriteFunctions.h>
#include <animaVelocityUtils.h>
#include <animaLinearTransformEstimationTools.h>

#include <animaAsymmetricBMRegistrationMethod.h>
#include <animaSymmetricBMRegistrationMethod.h>
#include <animaKissingSymmetricBMRegistrationMethod.h>
#include <animaResampleImageFilter.h>

#include <animaAnatomicalBlockMatcher.h>

namespace anima
{

template <unsigned int ImageDimension>
PyramidalDenseSVFMatchingBridge<ImageDimension>::PyramidalDenseSVFMatchingBridge()
{
    m_ReferenceImage = NULL;
    m_FloatingImage = NULL;

    m_OutputTransform = BaseTransformType::New();
    m_OutputTransform->SetIdentity();

    m_outputTransformFile = "";

    m_OutputImage = NULL;

    m_ReferenceMinimalValue = 0.0;
    m_FloatingMinimalValue = 0.0;

    m_BlockSize = 5;
    m_BlockSpacing = 2;
    m_StDevThreshold = 5;

    m_SymmetryType = Asymmetric;
    m_Transform = Translation;
    m_AffineDirection = 1;
    m_Metric = SquaredCorrelation;
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
    m_NumberOfPyramidLevels = 3;
    m_LastPyramidLevel = 0;
    m_PercentageKept = 0.8;
    this->SetNumberOfWorkUnits(itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads());

    m_Abort = false;
    m_Verbose = true;

    m_callback = itk::CStyleCommand::New();
    m_callback->SetClientData ((void *) this);
    m_callback->SetCallback (ManageProgress);
}

template <unsigned int ImageDimension>
PyramidalDenseSVFMatchingBridge<ImageDimension>::~PyramidalDenseSVFMatchingBridge()
{
}

template <unsigned int ImageDimension>
void
PyramidalDenseSVFMatchingBridge<ImageDimension>::Abort()
{
    m_Abort = true;

    if(m_bmreg)
        m_bmreg->Abort();
}

template <unsigned int ImageDimension>
void
PyramidalDenseSVFMatchingBridge<ImageDimension>::Update()
{
    m_Abort = false;

    // progress management
    m_progressReporter = new itk::ProgressReporter(this, 0, GetNumberOfPyramidLevels()*m_MaximumIterations);
    this->AddObserver(itk::ProgressEvent(), m_progressCallback);

    this->InvokeEvent(itk::StartEvent());

    // Compute minimal value of reference and floating images
    using MinMaxFilterType = itk::MinimumMaximumImageFilter <InputImageType>;
    typename MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();
    minMaxFilter->SetInput(m_ReferenceImage);
    if (this->GetNumberOfWorkUnits() != 0)
        minMaxFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    minMaxFilter->Update();

    m_ReferenceMinimalValue = minMaxFilter->GetMinimum();

    minMaxFilter = MinMaxFilterType::New();
    minMaxFilter->SetInput(m_FloatingImage);
    if (this->GetNumberOfWorkUnits() != 0)
        minMaxFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    minMaxFilter->Update();

    m_FloatingMinimalValue = minMaxFilter->GetMinimum();

    // Only CT images are below zero, little hack to set minimal values to either -1024 or 0
    if (m_ReferenceMinimalValue < 0.0)
        m_ReferenceMinimalValue = -1024;
    else
        m_ReferenceMinimalValue = 0.0;

    if (m_FloatingMinimalValue < 0.0)
        m_FloatingMinimalValue = -1024;
    else
        m_FloatingMinimalValue = 0.0;

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

        typename MaskImageType::Pointer maskGenerationImage = ITK_NULLPTR;
        if (m_BlockGenerationPyramid)
        {
            maskGenerationImage = m_BlockGenerationPyramid->GetOutput(i);
            maskGenerationImage->DisconnectPipeline();
        }

        // Update field to match the current resolution
        if (m_OutputTransform->GetParametersAsVectorField() != NULL)
        {
            typedef itk::ResampleImageFilter<VelocityFieldType,VelocityFieldType> VectorResampleFilterType;
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

        if (m_Verbose)
        {
            std::cout << "Processing pyramid level " << i << std::endl;
            std::cout << "Image size: " << refImage->GetLargestPossibleRegion().GetSize() << std::endl;
        }

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

            if (this->GetNumberOfWorkUnits() != 0)
                agreg->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

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

            if (this->GetNumberOfWorkUnits() != 0)
                agreg->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

            agreg->SetGeometryInformation(refImage.GetPointer());

            agregPtr = agreg;
        }

        agregPtr->SetVerboseAgregation(m_Verbose);

        // Init matcher
        typedef anima::AnatomicalBlockMatcher <InputImageType> BlockMatcherType;

        BlockMatcherType *mainMatcher = new BlockMatcherType;
        BlockMatcherType *reverseMatcher = 0;
        mainMatcher->SetBlockPercentageKept(GetPercentageKept());
        mainMatcher->SetBlockSize(GetBlockSize());
        mainMatcher->SetBlockSpacing(GetBlockSpacing());
        mainMatcher->SetBlockVarianceThreshold(GetStDevThreshold() * GetStDevThreshold());
        mainMatcher->SetBlockGenerationMask(maskGenerationImage);
        mainMatcher->SetDefaultBackgroundValue(m_FloatingMinimalValue);

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

                reverseMatcher = new BlockMatcherType;
                reverseMatcher->SetBlockPercentageKept(GetPercentageKept());
                reverseMatcher->SetBlockSize(GetBlockSize());
                reverseMatcher->SetBlockSpacing(GetBlockSpacing());
                reverseMatcher->SetBlockVarianceThreshold(GetStDevThreshold() * GetStDevThreshold());
                reverseMatcher->SetBlockGenerationMask(maskGenerationImage);
                reverseMatcher->SetDefaultBackgroundValue(m_ReferenceMinimalValue);
                reverseMatcher->SetVerbose(m_Verbose);

                tmpReg->SetReverseBlockMatcher(reverseMatcher);
                m_bmreg = tmpReg;
                break;
            }

            case Kissing:
            {
                typedef typename anima::KissingSymmetricBMRegistrationMethod <InputImageType> BlockMatchRegistrationType;
                typename BlockMatchRegistrationType::Pointer tmpReg = BlockMatchRegistrationType::New();
                tmpReg->SetReferenceBackgroundValue(m_ReferenceMinimalValue);
                tmpReg->SetFloatingBackgroundValue(m_FloatingMinimalValue);

                m_bmreg = tmpReg;
                break;
            }
        }

        mainMatcher->SetVerbose(m_Verbose);
        m_bmreg->SetBlockMatcher(mainMatcher);
        m_bmreg->SetBCHCompositionOrder(m_BCHCompositionOrder);
        m_bmreg->SetExponentiationOrder(m_ExponentiationOrder);

        if (m_progressCallback)
        {
            // we cannot connect directly bmreg to m_progressCallback
            // we need to create a new progressReporter with more iterations (m_progressReporter),
            // to listen to progress events from bmreg and to send new ones to m_progressCallback
            m_bmreg->AddObserver(itk::ProgressEvent(), m_callback);
        }

        m_bmreg->SetAgregator(agregPtr);

        if (this->GetNumberOfWorkUnits() != 0)
            m_bmreg->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

        m_bmreg->SetFixedImage(refImage);
        m_bmreg->SetMovingImage(floImage);

        m_bmreg->SetSVFElasticRegSigma(m_ElasticSigma * meanSpacing);

        typedef anima::ResampleImageFilter<InputImageType, InputImageType,
                                         typename BaseAgregatorType::ScalarType> ResampleFilterType;

        typename ResampleFilterType::Pointer refResampler = ResampleFilterType::New();
        refResampler->SetSize(floImage->GetLargestPossibleRegion().GetSize());
        refResampler->SetOutputOrigin(floImage->GetOrigin());
        refResampler->SetOutputSpacing(floImage->GetSpacing());
        refResampler->SetOutputDirection(floImage->GetDirection());
        refResampler->SetDefaultPixelValue(m_ReferenceMinimalValue);
        refResampler->SetNumberOfWorkUnits(GetNumberOfWorkUnits());
        m_bmreg->SetReferenceImageResampler(refResampler);

        typename ResampleFilterType::Pointer movingResampler = ResampleFilterType::New();
        movingResampler->SetSize(refImage->GetLargestPossibleRegion().GetSize());
        movingResampler->SetOutputOrigin(refImage->GetOrigin());
        movingResampler->SetOutputSpacing(refImage->GetSpacing());
        movingResampler->SetOutputDirection(refImage->GetDirection());
        movingResampler->SetDefaultPixelValue(m_FloatingMinimalValue);
        movingResampler->SetNumberOfWorkUnits(GetNumberOfWorkUnits());
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

            case Directional_Affine:
                mainMatcher->SetBlockTransformType(BlockMatcherType::Superclass::Directional_Affine);
                mainMatcher->SetAffineDirection(m_AffineDirection);
                if (reverseMatcher)
                {
                    reverseMatcher->SetBlockTransformType(BlockMatcherType::Superclass::Directional_Affine);
                    reverseMatcher->SetAffineDirection(m_AffineDirection);
                }
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

        switch (GetMetric())
        {
            case SquaredCorrelation:
                mainMatcher->SetSimilarityType(BlockMatcherType::SquaredCorrelation);
                if (reverseMatcher)
                    reverseMatcher->SetSimilarityType(BlockMatcherType::SquaredCorrelation);
                break;
            case Correlation:
                mainMatcher->SetSimilarityType(BlockMatcherType::Correlation);
                if (reverseMatcher)
                    reverseMatcher->SetSimilarityType(BlockMatcherType::Correlation);
                break;
            case MeanSquares:
            default:
                mainMatcher->SetSimilarityType(BlockMatcherType::MeanSquares);
                if (reverseMatcher)
                    reverseMatcher->SetSimilarityType(BlockMatcherType::MeanSquares);
                break;
        }

        m_bmreg->SetMaximumIterations(m_MaximumIterations);
        m_bmreg->SetMinimalTransformError(m_MinimalTransformError);
        m_bmreg->SetInitialTransform(m_OutputTransform.GetPointer());

        mainMatcher->SetNumberOfWorkUnits(GetNumberOfWorkUnits());
        mainMatcher->SetOptimizerMaximumIterations(GetOptimizerMaximumIterations());

        double sr = GetSearchRadius();
        mainMatcher->SetSearchRadius(sr);

        double sar = GetSearchAngleRadius();
        mainMatcher->SetSearchAngleRadius(sar);

        double scr = GetSearchScaleRadius();
        mainMatcher->SetSearchScaleRadius(scr);

        double fr = GetFinalRadius();
        mainMatcher->SetFinalRadius(fr);

        double ss = GetStepSize();
        mainMatcher->SetStepSize(ss);

        double tub = GetTranslateUpperBound();
        mainMatcher->SetTranslateMax(tub);

        double aub = GetAngleUpperBound();
        mainMatcher->SetAngleMax(aub);

        double scub = GetScaleUpperBound();
        mainMatcher->SetScaleMax(scub);

        if (reverseMatcher)
        {
            reverseMatcher->SetNumberOfWorkUnits(GetNumberOfWorkUnits());
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

        m_bmreg->SetVerboseProgression(m_Verbose);

        try
        {
            m_bmreg->Update();
        }
        catch( itk::ExceptionObject & err )
        {
            std::cout << "ExceptionObject caught !" << err << std::endl;
            exit(-1);
        }

        const BaseTransformType *resTrsf = dynamic_cast <const BaseTransformType *> (m_bmreg->GetOutput()->Get());
        m_OutputTransform->SetParametersAsVectorField(resTrsf->GetParametersAsVectorField());

        delete mainMatcher;
        if (reverseMatcher)
            delete reverseMatcher;
        if (agregPtr)
            delete agregPtr;
    }

    if (m_Abort)
        std::cout << "Process aborted" << std::endl;

    this->InvokeEvent(itk::EndEvent());

    if (m_SymmetryType == Kissing)
    {
        VelocityFieldType *finalTrsfField = const_cast <VelocityFieldType *> (m_OutputTransform->GetParametersAsVectorField());
        typedef itk::MultiplyImageFilter <VelocityFieldType,itk::Image <float,ImageDimension>, VelocityFieldType> MultiplyFilterType;

        typename MultiplyFilterType::Pointer fieldMultiplier = MultiplyFilterType::New();
        fieldMultiplier->SetInput(finalTrsfField);
        fieldMultiplier->SetConstant(2.0);
        fieldMultiplier->SetNumberOfWorkUnits(GetNumberOfWorkUnits());
        fieldMultiplier->InPlaceOn();

        fieldMultiplier->Update();

        VelocityFieldType *outputField = fieldMultiplier->GetOutput();
        m_OutputTransform->SetParametersAsVectorField(fieldMultiplier->GetOutput());
        outputField->DisconnectPipeline();
    }

    DisplacementFieldTransformPointer outputDispTrsf = DisplacementFieldTransformType::New();
    anima::GetSVFExponential(m_OutputTransform.GetPointer(), outputDispTrsf.GetPointer(), m_ExponentiationOrder, GetNumberOfWorkUnits(), false);

    typedef typename anima::ResampleImageFilter<InputImageType, InputImageType,
                                                typename BaseAgregatorType::ScalarType> ResampleFilterType;
    typename ResampleFilterType::Pointer tmpResample = ResampleFilterType::New();
    tmpResample->SetTransform(outputDispTrsf);
    tmpResample->SetInput(m_FloatingImage);

    tmpResample->SetSize(m_ReferenceImage->GetLargestPossibleRegion().GetSize());
    tmpResample->SetOutputOrigin(m_ReferenceImage->GetOrigin());
    tmpResample->SetOutputSpacing(m_ReferenceImage->GetSpacing());
    tmpResample->SetOutputDirection(m_ReferenceImage->GetDirection());
    tmpResample->SetDefaultPixelValue(m_ReferenceMinimalValue);
    tmpResample->Update();

    m_OutputImage = tmpResample->GetOutput();
    m_OutputImage->DisconnectPipeline();
}

template <unsigned int ImageDimension>
typename PyramidalDenseSVFMatchingBridge<ImageDimension>::DisplacementFieldTransformPointer
PyramidalDenseSVFMatchingBridge<ImageDimension>::GetOutputDisplacementFieldTransform()
{
    DisplacementFieldTransformPointer outputDispTrsf = DisplacementFieldTransformType::New();

    anima::GetSVFExponential(m_OutputTransform.GetPointer(), outputDispTrsf.GetPointer(), m_ExponentiationOrder, this->GetNumberOfWorkUnits(), false);

    return outputDispTrsf;
}

template <unsigned int ImageDimension>
void
PyramidalDenseSVFMatchingBridge<ImageDimension>::EmitProgress(int prog)
{
    if (m_progressReporter)
        m_progressReporter->CompletedPixel();
}

template <unsigned int ImageDimension>
void PyramidalDenseSVFMatchingBridge<ImageDimension>::ManageProgress (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    PyramidalDenseSVFMatchingBridge * source = reinterpret_cast<PyramidalDenseSVFMatchingBridge *> (clientData);
    itk::ProcessObject *processObject = (itk::ProcessObject *) caller;

    if (source && processObject)
        source->EmitProgress(processObject->GetProgress() * 100);
}

template <unsigned int ImageDimension>
void
PyramidalDenseSVFMatchingBridge<ImageDimension>::WriteOutputs()
{
    std::cout << "Writing output image to: " << m_resultFile << std::endl;
    anima::writeImage <InputImageType> (m_resultFile,m_OutputImage);

    if (m_outputTransformFile != "")
    {
        std::cout << "Writing output SVF to: " << m_outputTransformFile << std::endl;
        anima::writeImage <VelocityFieldType> (m_outputTransformFile,
                                               const_cast <VelocityFieldType *> (m_OutputTransform->GetParametersAsVectorField()));
    }
}

template <unsigned int ImageDimension>
void
PyramidalDenseSVFMatchingBridge<ImageDimension>::SetupPyramids()
{
    // Create pyramid here, check images actually are of the same size.
    m_ReferencePyramid = PyramidType::New();

    m_ReferencePyramid->SetInput(m_ReferenceImage);
    m_ReferencePyramid->SetNumberOfLevels(m_NumberOfPyramidLevels);

    if (this->GetNumberOfWorkUnits() != 0)
        m_ReferencePyramid->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

    typedef anima::ResampleImageFilter<InputImageType, InputImageType,
                                     typename BaseAgregatorType::ScalarType> ResampleFilterType;

    typename ResampleFilterType::Pointer refResampler = ResampleFilterType::New();
    refResampler->SetDefaultPixelValue(m_ReferenceMinimalValue);
    m_ReferencePyramid->SetImageResampler(refResampler);

    m_ReferencePyramid->Update();

    // Create pyramid for floating image
    m_FloatingPyramid = PyramidType::New();

    m_FloatingPyramid->SetInput(m_FloatingImage);
    m_FloatingPyramid->SetNumberOfLevels(m_NumberOfPyramidLevels);

    if (this->GetNumberOfWorkUnits() != 0)
        m_FloatingPyramid->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

    typename ResampleFilterType::Pointer floResampler = ResampleFilterType::New();
    floResampler->SetDefaultPixelValue(m_FloatingMinimalValue);
    m_FloatingPyramid->SetImageResampler(floResampler);

    m_FloatingPyramid->Update();

    m_BlockGenerationPyramid = 0;
    if (m_BlockGenerationMask)
    {
        typedef anima::ResampleImageFilter<MaskImageType, MaskImageType,
                typename BaseAgregatorType::ScalarType> MaskResampleFilterType;

        typename MaskResampleFilterType::Pointer maskResampler = MaskResampleFilterType::New();

        m_BlockGenerationPyramid = MaskPyramidType::New();
        m_BlockGenerationPyramid->SetImageResampler(maskResampler);
        m_BlockGenerationPyramid->SetInput(m_BlockGenerationMask);
        m_BlockGenerationPyramid->SetNumberOfLevels(GetNumberOfPyramidLevels());
        m_BlockGenerationPyramid->SetNumberOfWorkUnits(GetNumberOfWorkUnits());
        m_BlockGenerationPyramid->Update();
    }
}

} // end of namespace
