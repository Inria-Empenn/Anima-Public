#pragma once

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkMultiThreader.h>
#include <itkCenteredTransformInitializer.h>

#include <animaAsymmetricBlockMatchingRegistrationMethod.h>
#include <animaSymmetricBlockMatchingRegistrationMethod.h>
#include <animaKissingSymmetricBlockMatchingRegistrationMethod.h>

#include <animaLSWTransformAgregator.h>
#include <animaLTSWTransformAgregator.h>
#include <animaMEstTransformAgregator.h>

// ------------------------------

namespace anima
{

template <unsigned int ImageDimension>
PyramidalBlockMatchingBridge<ImageDimension>::PyramidalBlockMatchingBridge()
{
    m_InitialTransform = NULL;
    m_ReferenceImage = NULL;
    m_FloatingImage = NULL;

    m_OutputTransform = NULL;
    m_outputTransformFile = "";

    m_OutputImage = NULL;

    m_BlockSize = 5;
    m_BlockSpacing = 5;
    m_StDevThreshold = 5;

    m_SymmetryType = Asymmetric;
    m_Transform = Translation;
    m_Metric = SquaredCorrelation;
    m_Optimizer = Bobyqa;

    m_MaximumIterations = 10;
    m_MinimalTransformError = 0.01;
    m_OptimizerMaximumIterations = 100;
    m_SearchRadius = 2;
    m_SearchAngleRadius = 5;
    m_SearchSkewRadius = 5;
    m_SearchScaleRadius = 0.1;
    m_FinalRadius = 0.001;
    m_StepSize = 1;
    m_TranslateUpperBound = 50;
    m_AngleUpperBound = 180;
    m_SkewUpperBound = 45;
    m_ScaleUpperBound = 3;
    m_Agregator = MEstimation;
    m_OutputTransformType = outRigid;
    m_AgregThreshold = 0.5;
    m_SeStoppingThreshold = 0.01;
    m_NumberOfPyramidLevels = 3;
    m_LastPyramidLevel = 0;
    m_PercentageKept = 0.8;
    m_InitializeOnCenterOfGravity = true;

    this->SetNumberOfThreads(itk::MultiThreader::GetGlobalDefaultNumberOfThreads());

    m_Abort = false;

    m_callback = itk::CStyleCommand::New();
    m_callback->SetClientData ((void *) this);
    m_callback->SetCallback (ManageProgress);
}

template <unsigned int ImageDimension>
PyramidalBlockMatchingBridge<ImageDimension>::~PyramidalBlockMatchingBridge()
{
}

template <unsigned int ImageDimension>
void PyramidalBlockMatchingBridge<ImageDimension>::Abort()
{
    m_Abort = true;

    if(m_bmreg)
        m_bmreg->Abort();
}

template <unsigned int ImageDimension>
void PyramidalBlockMatchingBridge<ImageDimension>::SetInitialTransform (std::string initialTransformFile)
{
    if (initialTransformFile != "")
    {
        itk::TransformFileReader::Pointer tmpTrRead = itk::TransformFileReader::New();
        tmpTrRead->SetFileName (initialTransformFile);

        try
        {
            tmpTrRead->Update();

            itk::TransformFileReader::TransformListType trsfList = * (tmpTrRead->GetTransformList());
            itk::TransformFileReader::TransformListType::iterator tr_it = trsfList.begin();

            m_InitialTransform = dynamic_cast <AffineTransformType *> ((*tr_it).GetPointer());

        }
        catch (itk::ExceptionObject &e)
        {
            std::cerr << "Unable to read initial transform... Exiting..." << std::endl;
        }
    }
}

template <unsigned int ImageDimension>
void PyramidalBlockMatchingBridge<ImageDimension>::InitializeBlocksOnImage(InitializerPointer &initPtr, InputImageType *image)
{
    // Init blocks
    initPtr = InitializerType::New();

    initPtr->AddReferenceImage(image);
    initPtr->SetNumberOfThreads(GetNumberOfThreads());

    initPtr->SetPercentageKept(GetPercentageKept());
    initPtr->SetBlockSize(GetBlockSize());
    initPtr->SetBlockSpacing(GetBlockSpacing());
    initPtr->SetScalarVarianceThreshold(GetStDevThreshold() * GetStDevThreshold());

    initPtr->SetRequestedRegion(image->GetLargestPossibleRegion());
}

template <unsigned int ImageDimension>
void PyramidalBlockMatchingBridge<ImageDimension>::Update()
{
    typedef BaseTransformAgregator<ImageDimension> BaseAgreg;

    m_Abort = false;

    // progress management
    m_progressReporter = new itk::ProgressReporter(this, 0, GetNumberOfPyramidLevels()*this->m_MaximumIterations);
    this->AddObserver(itk::ProgressEvent(), m_progressCallback);

    this->InvokeEvent(itk::StartEvent());

    this->SetupPyramids();

    // Iterate over pyramid levels
    for (unsigned int i = 0;i < GetNumberOfPyramidLevels() && !m_Abort; ++i)
    {
        if (i + GetLastPyramidLevel() >= m_ReferencePyramid->GetNumberOfLevels())
            continue;

        typename InputImageType::Pointer refImage = m_ReferencePyramid->GetOutput(i);
        refImage->DisconnectPipeline();

        typename InputImageType::Pointer floImage = m_FloatingPyramid->GetOutput(i);
        floImage->DisconnectPipeline();

        std::cout << "Processing pyramid level " << i << std::endl;
        std::cout << "Image size: " << refImage->GetLargestPossibleRegion().GetSize() << std::endl;

        // Init matcher
        switch (m_SymmetryType)
        {
            case Asymmetric:
            {
                typedef typename anima::AsymmetricBlockMatchingRegistrationMethod <InputImageType> BlockMatchRegistrationType;
                m_bmreg = BlockMatchRegistrationType::New();

                typename InitializerType::Pointer initPtr;
                this->InitializeBlocksOnImage(initPtr, refImage);

                m_bmreg->SetBlockRegions(initPtr->GetOutput());
                std::cout << "Generated " << initPtr->GetOutput().size() << " blocks..." << std::endl;

                break;
            }

            case Symmetric:
            {
                typedef typename anima::SymmetricBlockMatchingRegistrationMethod <InputImageType> BlockMatchRegistrationType;
                typename BlockMatchRegistrationType::Pointer tmpReg = BlockMatchRegistrationType::New();

                typename InitializerType::Pointer initPtr;
                this->InitializeBlocksOnImage(initPtr, refImage);

                tmpReg->SetFixedBlockRegions(initPtr->GetOutput());
                std::cout << "Generated " << initPtr->GetOutput().size() << " blocks..." << std::endl;

                this->InitializeBlocksOnImage(initPtr, floImage);

                tmpReg->SetMovingBlockRegions(initPtr->GetOutput());
                std::cout << "Generated " << initPtr->GetOutput().size() << " blocks..." << std::endl;

                m_bmreg = tmpReg;
                break;
            }

            case Kissing:
            {
                typedef typename anima::KissingSymmetricBlockMatchingRegistrationMethod <InputImageType> BlockMatchRegistrationType;
                typename BlockMatchRegistrationType::Pointer tmpReg = BlockMatchRegistrationType::New();

                tmpReg->SetBlockPercentageKept(GetPercentageKept());
                tmpReg->SetBlockSize(GetBlockSize());
                tmpReg->SetBlockSpacing(GetBlockSpacing());

                m_bmreg = tmpReg;

                break;
            }
        }

        if (m_progressCallback)
        {
            // we cannot connect directly bmreg to m_progressCallback
            // we need to create a new progressReporter with more iterations (m_progressReporter),
            // to listen to progress events from bmreg and to send new ones to m_progressCallback
            m_bmreg->AddObserver(itk::ProgressEvent(), m_callback);
        }

        m_bmreg->SetNumberOfThreads(GetNumberOfThreads());
        m_bmreg->SetBlockScalarVarianceThreshold(GetStDevThreshold() * GetStDevThreshold());

        m_bmreg->SetFixedImage(refImage);
        m_bmreg->SetMovingImage(floImage);

        BaseAgreg* agreg = NULL;

        switch (GetAgregator())
        {
            case LeastSquares:
            {
                typedef LSWTransformAgregator<ImageDimension> Agreg;
                agreg = new Agreg;
                break;
            }

            case LeastTrimmedSquares:
            {
                typedef LTSWTransformAgregator<ImageDimension> Agreg;
                Agreg *tmpAg = new Agreg;

                tmpAg->SetLTSCut(GetAgregThreshold());
                tmpAg->SeStoppingThreshold(GetSeStoppingThreshold());

                agreg = tmpAg;
                break;
            }

            case MEstimation:
            {
                typedef MEstTransformAgregator<ImageDimension> Agreg;
                Agreg *tmpAg = new Agreg;

                tmpAg->SetMEstimateFactor(GetAgregThreshold());
                tmpAg->SeStoppingThreshold(GetSeStoppingThreshold());

                agreg = tmpAg;
                break;
            }
        }

        switch (GetOutputTransformType())
        {
            case outTranslation:
                agreg->SetOutputTransformType(BaseAgreg::TRANSLATION);
                break;
            case outRigid:
                agreg->SetOutputTransformType(BaseAgreg::RIGID);
                break;
            case outAffine:
            default:
                agreg->SetOutputTransformType(BaseAgreg::AFFINE);
                break;
        }

        m_bmreg->SetAgregator(agreg);

        switch (GetTransform())
        {
            case Translation:
                m_bmreg->SetTransformKind(BaseBlockMatchRegistrationType::Translation);
                break;
            case Rigid:
                m_bmreg->SetTransformKind(BaseBlockMatchRegistrationType::Rigid);
                break;
            case Affine:
            default:
                m_bmreg->SetTransformKind(BaseBlockMatchRegistrationType::Affine);
                break;
        }

        switch (GetOptimizer())
        {
            case Exhaustive:
                m_bmreg->SetOptimizerKind(BaseBlockMatchRegistrationType::Exhaustive);
                break;
            case Bobyqa:
            default:
                m_bmreg->SetOptimizerKind(BaseBlockMatchRegistrationType::Bobyqa);
                break;
        }

        switch (GetMetric())
        {
            case Correlation:
                m_bmreg->SetMetricKind(BaseBlockMatchRegistrationType::Correlation);
                break;
            case SquaredCorrelation:
                m_bmreg->SetMetricKind(BaseBlockMatchRegistrationType::SquaredCorrelation);
                break;
            case MeanSquares:
            default:
                m_bmreg->SetMetricKind(BaseBlockMatchRegistrationType::MeanSquares);
                break;
        }

        m_bmreg->SetMaximumIterations(GetMaximumIterations());
        m_bmreg->SetOptimizerMaximumIterations(GetOptimizerMaximumIterations());
        m_bmreg->SetMinimalTransformError(GetMinimalTransformError());

        m_bmreg->SetInitialTransform(m_OutputTransform);

        double sr = GetSearchRadius();
        m_bmreg->SetSearchRadius(sr);

        double sar = GetSearchAngleRadius();
        m_bmreg->SetSearchAngleRadius(sar);

        double skr = GetSearchSkewRadius();
        m_bmreg->SetSearchSkewRadius(skr);

        double scr = GetSearchScaleRadius();
        m_bmreg->SetSearchScaleRadius(scr);

        double fr = GetFinalRadius();
        m_bmreg->SetFinalRadius(fr);

        double ss = GetStepSize();
        m_bmreg->SetStepSize(ss);

        double tub = GetTranslateUpperBound();
        m_bmreg->SetTranslateMax(tub);

        double aub = GetAngleUpperBound();
        m_bmreg->SetAngleMax(aub);

        double skub = GetSkewUpperBound();
        m_bmreg->SetSkewMax(skub);

        double scub = GetScaleUpperBound();
        m_bmreg->SetScaleMax(scub);

        try
        {
            m_bmreg->Update();
            std::cout << "Block Matching Registration stop condition "
                      << m_bmreg->GetStopConditionDescription()
                      << std::endl;
        }
        catch( itk::ExceptionObject & err )
        {
            std::cout << "ExceptionObject caught in bmreg startregistration ! " << err << std::endl;
            exit(-1);
        }

        // Polyrigid will have to be handled here
        AffineTransformType *tmpTrsf = dynamic_cast<AffineTransformType *>(m_OutputTransform.GetPointer());
        tmpTrsf->SetParameters(m_bmreg->GetOutput()->Get()->GetParameters());
    }

    if (m_Abort)
        std::cout << "Process aborted" << std::endl;

    this->InvokeEvent(itk::EndEvent());

    AffineTransformType *tmpTrsf = dynamic_cast<AffineTransformType *>(m_OutputTransform.GetPointer());
    
    if (m_SymmetryType == Kissing)
    {
        typename AffineTransformType::Pointer trsfCopy = AffineTransformType::New();
        trsfCopy->SetMatrix(tmpTrsf->GetMatrix());
        trsfCopy->SetOffset(tmpTrsf->GetOffset());

        tmpTrsf->Compose(trsfCopy);
    }

    if (!m_InitialTransform.IsNull())
        tmpTrsf->Compose(m_InitialTransform, false);

    typedef typename anima::ResampleImageFilter<InputImageType, InputImageType,
                                                typename AgregatorType::ScalarType> ResampleFilterType;
    typename ResampleFilterType::Pointer tmpResample = ResampleFilterType::New();
    tmpResample->SetTransform(m_OutputTransform);
    tmpResample->SetInput(m_FloatingImage);

    tmpResample->SetSize(m_ReferenceImage->GetLargestPossibleRegion().GetSize());
    tmpResample->SetOutputOrigin(m_ReferenceImage->GetOrigin());
    tmpResample->SetOutputSpacing(m_ReferenceImage->GetSpacing());
    tmpResample->SetOutputDirection(m_ReferenceImage->GetDirection());
    tmpResample->SetDefaultPixelValue(0);
    tmpResample->Update();

    m_OutputImage = tmpResample->GetOutput();
}

template <unsigned int ImageDimension>
void PyramidalBlockMatchingBridge<ImageDimension>::EmitProgress(int prog)
{
   if (m_progressReporter)
       m_progressReporter->CompletedPixel();
}

template <unsigned int ImageDimension>
void PyramidalBlockMatchingBridge<ImageDimension>::ManageProgress (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    PyramidalBlockMatchingBridge * source = reinterpret_cast<PyramidalBlockMatchingBridge *> (clientData);
    itk::ProcessObject *processObject = (itk::ProcessObject *) caller;

    if (source && processObject)
        source->EmitProgress(processObject->GetProgress() * 100);
}

template <unsigned int ImageDimension>
void PyramidalBlockMatchingBridge<ImageDimension>::WriteOutputs()
{
    std::cout << "Writing output image to: " << GetResultFile() << std::endl;

    typename itk::ImageFileWriter <InputImageType>::Pointer imageWriter = itk::ImageFileWriter <InputImageType>::New();
    imageWriter->SetUseCompression(true);
    imageWriter->SetInput(m_OutputImage);
    imageWriter->SetFileName(GetResultFile());

    imageWriter->Update();

    if (GetOutputTransformFile() != "")
    {
        std::cout << "Writing output transform to: " << GetOutputTransformFile() << std::endl;
        itk::TransformFileWriter::Pointer writer = itk::TransformFileWriter::New();
        writer->SetInput(m_OutputTransform);
        writer->SetFileName(GetOutputTransformFile());
        writer->Update();
    }
}

template <unsigned int ImageDimension>
        void PyramidalBlockMatchingBridge<ImageDimension>::SetupPyramids()
{
    // Create pyramid here, check images actually are of the same size.
    typedef anima::ResampleImageFilter<InputImageType, InputImageType,
                                     typename AgregatorType::ScalarType> ResampleFilterType;
    typedef typename itk::CenteredTransformInitializer<AffineTransformType, InputImageType, InputImageType> TransformInitializerType;

    m_ReferencePyramid = PyramidType::New();

    m_ReferencePyramid->SetInput(m_ReferenceImage);
    m_ReferencePyramid->SetNumberOfLevels(GetNumberOfPyramidLevels());
    m_ReferencePyramid->SetNumberOfThreads(GetNumberOfThreads());
    m_ReferencePyramid->Update();

    InputImagePointer initialFloatingImage = const_cast <InputImageType *> (m_FloatingImage.GetPointer());

    // Compute initial transform if needed to get a decent initial floating image
    if (!m_InitialTransform.IsNull())
    {
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
        initialFloatingImage->DisconnectPipeline();
    }
    else
    {
        m_InitialTransform = NULL;

        m_InitialTransform = AffineTransformType::New();
        m_InitialTransform->SetIdentity();

        if (m_InitializeOnCenterOfGravity)
        {
            typename TransformInitializerType::Pointer initializer = TransformInitializerType::New();
            initializer->SetTransform(m_InitialTransform);
            initializer->SetFixedImage(m_ReferenceImage);
            initializer->SetMovingImage(m_FloatingImage);
            initializer->MomentsOn();
            initializer->InitializeTransform();
        }

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
        initialFloatingImage->DisconnectPipeline();
    }

    // Create pyramid for floating image
    m_FloatingPyramid = PyramidType::New();

    m_FloatingPyramid->SetInput(initialFloatingImage);
    m_FloatingPyramid->SetNumberOfLevels(GetNumberOfPyramidLevels());
    m_FloatingPyramid->SetNumberOfThreads(GetNumberOfThreads());
    m_FloatingPyramid->Update();
}

} // end of namespace anima
