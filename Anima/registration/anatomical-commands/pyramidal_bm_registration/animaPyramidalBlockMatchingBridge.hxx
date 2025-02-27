#pragma once
#include "animaPyramidalBlockMatchingBridge.h"

#include <animaReadWriteFunctions.h>
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkCenteredTransformInitializer.h>
#include <itkMinimumMaximumImageFilter.h>
#include <itkAddImageFilter.h>

#include <animaAsymmetricBMRegistrationMethod.h>
#include <animaSymmetricBMRegistrationMethod.h>
#include <animaKissingSymmetricBMRegistrationMethod.h>

#include <animaAnatomicalBlockMatcher.h>

#include <animaLSWTransformAgregator.h>
#include <animaLTSWTransformAgregator.h>
#include <animaMEstTransformAgregator.h>

namespace anima
{

template <unsigned int ImageDimension>
PyramidalBlockMatchingBridge<ImageDimension>::PyramidalBlockMatchingBridge()
{
    m_InitialTransform = nullptr;
    m_DirectionTransform = nullptr;
    m_ReferenceImage = nullptr;
    m_FloatingImage = nullptr;

    m_OutputTransform = nullptr;
    m_outputTransformFile = "";
    m_outputNearestRigidTransformFile = "";
    m_outputNearestSimilarityTransformFile = "";

    m_OutputImage = nullptr;

    m_ReferenceMinimalValue = 0.0;
    m_FloatingMinimalValue = 0.0;
    m_RegistrationPointLocation = 0.5;

    m_BlockSize = 5;
    m_BlockSpacing = 5;
    m_StDevThreshold = 5;

    m_SymmetryType = Asymmetric;
    m_Transform = Translation;
    m_AffineDirection = 1;
    m_Metric = SquaredCorrelation;
    m_Optimizer = Bobyqa;

    m_MaximumIterations = 10;
    m_MinimalTransformError = 0.01;
    m_OptimizerMaximumIterations = 100;
    m_StepSize = 1;
    m_TranslateUpperBound = 50;
    m_AngleUpperBound = 180;
    m_ScaleUpperBound = 3;
    m_Agregator = MEstimation;
    m_OutputTransformType = outRigid;
    m_AgregThreshold = 0.5;
    m_SeStoppingThreshold = 0.01;
    m_NumberOfPyramidLevels = 3;
    m_LastPyramidLevel = 0;
    m_PercentageKept = 0.8;
    m_TransformInitializationType = ClosestTransform;

    this->SetNumberOfWorkUnits(itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads());

    m_Abort = false;
    m_Verbose = true;

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
void PyramidalBlockMatchingBridge<ImageDimension>::SetDirectionTransform(std::string directionTransformFile)
{
    if (directionTransformFile != "")
    {
        itk::TransformFileReader::Pointer tmpTrRead = itk::TransformFileReader::New();
        tmpTrRead->SetFileName(directionTransformFile);

        try
        {
            tmpTrRead->Update();

            itk::TransformFileReader::TransformListType trsfList = *(tmpTrRead->GetTransformList());
            itk::TransformFileReader::TransformListType::iterator tr_it = trsfList.begin();

            m_DirectionTransform = dynamic_cast <AffineTransformType *> ((*tr_it).GetPointer());

        }
        catch (itk::ExceptionObject &e)
        {
            std::cerr << "Unable to read direction transform... Exiting..." << std::endl;
        }
    }
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

    // Compute minimal value of reference and Floating images
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

    // Set up pyramids of images and masks
    this->SetupPyramids();

    bool invertInputs = (m_RegistrationPointLocation < 0.5) && (m_SymmetryType == Kissing);
    if (invertInputs)
    {
        m_RegistrationPointLocation = 1.0 - m_RegistrationPointLocation;
        double floValue = m_FloatingMinimalValue;
        m_FloatingMinimalValue = m_ReferenceMinimalValue;
        m_ReferenceMinimalValue = floValue;
    }

    typedef anima::AnatomicalBlockMatcher <InputImageType> BlockMatcherType;

    // Iterate over pyramid levels
    for (unsigned int i = 0;i < GetNumberOfPyramidLevels() && !m_Abort; ++i)
    {
        if (i + GetLastPyramidLevel() >= m_ReferencePyramid->GetNumberOfLevels())
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

        BlockMatcherType *mainMatcher = new BlockMatcherType;
        BlockMatcherType *reverseMatcher = 0;
        mainMatcher->SetBlockPercentageKept(GetPercentageKept());
        mainMatcher->SetBlockSize(GetBlockSize());
        mainMatcher->SetBlockSpacing(GetBlockSpacing());
        mainMatcher->SetBlockVarianceThreshold(GetStDevThreshold() * GetStDevThreshold());
        mainMatcher->SetBlockGenerationMask(maskGenerationImage);
        mainMatcher->SetDefaultBackgroundValue(m_FloatingMinimalValue);

        if (m_Verbose)
        {
            std::cout << "Processing pyramid level " << i << std::endl;
            std::cout << "Image size: " << refImage->GetLargestPossibleRegion().GetSize() << std::endl;
        }

        // Init bm registration method
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
                reverseMatcher->SetVerbose(m_Verbose);
                reverseMatcher->SetBlockGenerationMask(maskGenerationImage);
                reverseMatcher->SetDefaultBackgroundValue(m_ReferenceMinimalValue);

                tmpReg->SetReverseBlockMatcher(reverseMatcher);
                m_bmreg = tmpReg;
                break;
            }

            case Kissing:
            default:
            {
                typedef typename anima::KissingSymmetricBMRegistrationMethod <InputImageType> BlockMatchRegistrationType;
                typename BlockMatchRegistrationType::Pointer tmpReg = BlockMatchRegistrationType::New();
                tmpReg->SetReferenceBackgroundValue(m_ReferenceMinimalValue);
                tmpReg->SetFloatingBackgroundValue(m_FloatingMinimalValue);
                tmpReg->SetRegistrationPointLocation(m_RegistrationPointLocation);

                m_bmreg = tmpReg;
                break;
            }
        }

        mainMatcher->SetVerbose(m_Verbose);
        m_bmreg->SetBlockMatcher(mainMatcher);

        if (m_progressCallback)
        {
            // we cannot connect directly bmreg to m_progressCallback
            // we need to create a new progressReporter with more iterations (m_progressReporter),
            // to listen to progress events from bmreg and to send new ones to m_progressCallback
            m_bmreg->AddObserver(itk::ProgressEvent(), m_callback);
        }

        if (this->GetNumberOfWorkUnits() != 0)
            m_bmreg->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

        m_bmreg->SetFixedImage(refImage);
        m_bmreg->SetMovingImage(floImage);

        typedef anima::ResampleImageFilter<InputImageType, InputImageType,
                typename AgregatorType::ScalarType> ResampleFilterType;

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

        BaseAgreg *agreg = NULL;

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
            case outAnisotropic_Sim:
                agreg->SetOutputTransformType(BaseAgreg::ANISOTROPIC_SIM);
                if (m_DirectionTransform)
                    agreg->SetOrthogonalDirectionMatrix(m_DirectionTransform->GetMatrix());
                break;
            case outAffine:
            default:
                agreg->SetOutputTransformType(BaseAgreg::AFFINE);
                break;
        }

        agreg->SetVerboseAgregation(m_Verbose);
        m_bmreg->SetAgregator(agreg);

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


        switch(GetOptimizer())
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
            case Correlation:
                mainMatcher->SetSimilarityType(BlockMatcherType::Correlation);
                if (reverseMatcher)
                    reverseMatcher->SetSimilarityType(BlockMatcherType::Correlation);
                break;
            case SquaredCorrelation:
                mainMatcher->SetSimilarityType(BlockMatcherType::SquaredCorrelation);
                if (reverseMatcher)
                    reverseMatcher->SetSimilarityType(BlockMatcherType::SquaredCorrelation);
                break;
            case MeanSquares:
            default:
                mainMatcher->SetSimilarityType(BlockMatcherType::MeanSquares);
                if (reverseMatcher)
                    reverseMatcher->SetSimilarityType(BlockMatcherType::MeanSquares);
                break;
        }

        m_bmreg->SetMaximumIterations(GetMaximumIterations());
        m_bmreg->SetMinimalTransformError(GetMinimalTransformError());
        m_bmreg->SetInitialTransform(m_OutputTransform);

        mainMatcher->SetNumberOfWorkUnits(GetNumberOfWorkUnits());
        mainMatcher->SetOptimizerMaximumIterations(GetOptimizerMaximumIterations());

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
            std::cerr << "ExceptionObject caught in bmreg startregistration ! " << err << std::endl;
            exit(-1);
        }

        if ((GetOutputTransformType() == outAnisotropic_Sim)||(GetOutputTransformType() == outAffine))
            m_EstimationBarycenter = agreg->GetEstimationBarycenter();

        // Polyrigid will have to be handled here
        AffineTransformType *tmpTrsf = dynamic_cast<AffineTransformType *>(m_OutputTransform.GetPointer());
        tmpTrsf->SetParameters(m_bmreg->GetOutput()->Get()->GetParameters());

        delete mainMatcher;
        if (reverseMatcher)
            delete reverseMatcher;
        if (agreg)
            delete agreg;
    }

    if (m_Abort)
        std::cout << "Process aborted" << std::endl;

    this->InvokeEvent(itk::EndEvent());
    this->RemoveAllObservers();

    AffineTransformType *tmpTrsf = dynamic_cast<AffineTransformType *>(m_OutputTransform.GetPointer());
    
    if (m_SymmetryType == Kissing)
    {
        unsigned int NDimensions = InputImageType::ImageDimension;
        vnl_matrix <double> affMatrix(NDimensions+1,NDimensions+1,0);
        affMatrix.set_identity();

        for (unsigned int i = 0;i < NDimensions;++i)
        {
            for (unsigned int j = 0;j < NDimensions;++j)
                affMatrix(i,j) = tmpTrsf->GetMatrix()(i,j);

            affMatrix(i,NDimensions) = tmpTrsf->GetOffset()[i];
        }

        affMatrix = anima::GetLogarithm(affMatrix);

        double multiplier = 1.0 / m_RegistrationPointLocation;
        if (invertInputs)
            multiplier *= -1.0;

        affMatrix *= multiplier;
        affMatrix = anima::GetExponential(affMatrix);

        vnl_matrix <double> affResult(NDimensions,NDimensions,0);
        typename AffineTransformType::OffsetType offset(NDimensions);

        for (unsigned int i = 0;i < NDimensions;++i)
        {
            for (unsigned int j = 0;j < NDimensions;++j)
                affResult(i,j) = affMatrix(i,j);

            offset[i] = affMatrix(i,NDimensions);
        }

        tmpTrsf->SetMatrix(affResult);
        tmpTrsf->SetOffset(offset);
    }

    if (!m_InitialTransform.IsNull())
        tmpTrsf->Compose(m_InitialTransform, false);

    if (invertInputs)
    {
        m_RegistrationPointLocation = 1.0 - m_RegistrationPointLocation;
        double floValue = m_FloatingMinimalValue;
        m_FloatingMinimalValue = m_ReferenceMinimalValue;
        m_ReferenceMinimalValue = floValue;
    }

    typedef typename anima::ResampleImageFilter<InputImageType, InputImageType, typename AgregatorType::ScalarType> ResampleFilterType;
    typename ResampleFilterType::Pointer tmpResample = ResampleFilterType::New();
    tmpResample->SetTransform(m_OutputTransform);
    tmpResample->SetInput(m_FloatingImage);

    tmpResample->SetSize(m_ReferenceImage->GetLargestPossibleRegion().GetSize());
    tmpResample->SetOutputOrigin(m_ReferenceImage->GetOrigin());
    tmpResample->SetOutputSpacing(m_ReferenceImage->GetSpacing());
    tmpResample->SetOutputDirection(m_ReferenceImage->GetDirection());
    tmpResample->SetDefaultPixelValue(m_FloatingMinimalValue);
    tmpResample->Update();

    m_OutputImage = tmpResample->GetOutput();
    m_OutputImage->DisconnectPipeline();
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
    anima::writeImage <InputImageType> (GetResultFile(),m_OutputImage);

    if (GetOutputTransformFile() != "")
    {
        std::cout << "Writing output transform to: " << GetOutputTransformFile() << std::endl;
        itk::TransformFileWriter::Pointer writer = itk::TransformFileWriter::New();
        writer->SetInput(m_OutputTransform);
        writer->SetFileName(GetOutputTransformFile());
        writer->Update();
    }

    if ((m_outputNearestRigidTransformFile == "")&&(m_outputNearestSimilarityTransformFile == ""))
        return;

    AffineTransformType *tmpTrsf = dynamic_cast<AffineTransformType *>(m_OutputTransform.GetPointer());

    typename AffineTransformType::MatrixType linearMatrix = tmpTrsf->GetMatrix();
    typename AffineTransformType::OffsetType transformOffset = tmpTrsf->GetOffset();
    vnl_svd<typename AffineTransformType::MatrixType::ValueType> UWVLinearMatrixSVD(linearMatrix.GetVnlMatrix().as_matrix());
    vnl_matrix<typename AffineTransformType::MatrixType::ValueType> leftRot = UWVLinearMatrixSVD.U()*vnl_determinant(UWVLinearMatrixSVD.U());
    vnl_matrix<typename AffineTransformType::MatrixType::ValueType> rightRot = UWVLinearMatrixSVD.V()*vnl_determinant(UWVLinearMatrixSVD.V());
    vnl_diag_matrix<typename AffineTransformType::MatrixType::ValueType> scal = UWVLinearMatrixSVD.W();

    PointType ybar;
    double isoScal = 0;
    for (unsigned int i = 0;i < ImageDimension;++i)
    {
        isoScal += std::abs(scal(i, i)) / ImageDimension;
        ybar[i] = transformOffset[i];

        for (unsigned int j = 0;j < ImageDimension;++j)
            ybar[i] += linearMatrix(i, j) * m_EstimationBarycenter[j];
    }

    typename AffineTransformType::OffsetType rigidOffset;
    typename AffineTransformType::MatrixType linearPartRigid;
    typename AffineTransformType::OffsetType similarityOffset;
    typename AffineTransformType::MatrixType linearPartSimilarity;
    typename AffineTransformType::OffsetType UOffset;
    linearPartRigid.Fill(0);
    linearPartSimilarity.Fill(0);

    for (unsigned int i = 0;i < ImageDimension;++i)
    {
        rigidOffset[i] = ybar[i];
        similarityOffset[i] = ybar[i];
        UOffset[i] = 0;
        for (unsigned int j = 0;j < ImageDimension; ++j)
        {
            for (unsigned int k = 0;k < ImageDimension;++k)
            {
                linearPartRigid(i, j) += leftRot(i, k) * rightRot(j, k);
                linearPartSimilarity(i, j) += isoScal * leftRot(i, k) * rightRot(j, k);
            }

            rigidOffset[i] -= linearPartRigid(i, j) * m_EstimationBarycenter[j];
            similarityOffset[i] -= linearPartSimilarity(i, j) * m_EstimationBarycenter[j];
        }
    }

    if (m_outputNearestRigidTransformFile != "")
    {
        AffineTransformPointer rigidTransform = AffineTransformType::New();
        rigidTransform->SetMatrix(linearPartRigid);
        rigidTransform->SetOffset(rigidOffset);
        itk::TransformFileWriter::Pointer rigidWriter = itk::TransformFileWriter::New();
        rigidWriter->SetInput(rigidTransform);
        rigidWriter->SetFileName(GetOutputNearestRigidTransformFile());
        std::cout << "Writing output nearest rigid transform to: " << GetOutputNearestRigidTransformFile() << std::endl;
        rigidWriter->Update();
    }

    if (m_outputNearestSimilarityTransformFile != "")
    {
        AffineTransformPointer similarityTransform = AffineTransformType::New();
        similarityTransform->SetMatrix(linearPartSimilarity);
        similarityTransform->SetOffset(similarityOffset);
        itk::TransformFileWriter::Pointer similarityWriter = itk::TransformFileWriter::New();
        similarityWriter->SetInput(similarityTransform);
        similarityWriter->SetFileName(GetOutputNearestSimilarityTransformFile());
        std::cout << "Writing output nearest similarity transform to: " << GetOutputNearestSimilarityTransformFile() << std::endl;
        similarityWriter->Update();
    }
}

template <unsigned int ImageDimension>
void PyramidalBlockMatchingBridge<ImageDimension>::SetupPyramids()
{
    // Create pyramid here, check images actually are of the same size.
    typedef anima::ResampleImageFilter<InputImageType, InputImageType,
            typename AgregatorType::ScalarType> ResampleFilterType;

    InputImagePointer initialFloatingImage = const_cast <InputImageType *> (m_FloatingImage.GetPointer());

    // Compute initial transform if needed to get a decent initial Floating image
    if (m_InitialTransform.IsNotNull())
    {
        typename ResampleFilterType::Pointer tmpResample = ResampleFilterType::New();
        tmpResample->SetTransform(m_InitialTransform);
        tmpResample->SetInput(m_FloatingImage);

        tmpResample->SetSize(m_ReferenceImage->GetLargestPossibleRegion().GetSize());
        tmpResample->SetOutputOrigin(m_ReferenceImage->GetOrigin());
        tmpResample->SetOutputSpacing(m_ReferenceImage->GetSpacing());
        tmpResample->SetOutputDirection(m_ReferenceImage->GetDirection());
        tmpResample->SetDefaultPixelValue(m_FloatingMinimalValue);
        tmpResample->Update();

        initialFloatingImage = tmpResample->GetOutput();
        initialFloatingImage->DisconnectPipeline();
    }
    else
    {
        m_InitialTransform = NULL;

        m_InitialTransform = AffineTransformType::New();
        m_InitialTransform->SetIdentity();

        typedef itk::ImageMomentsCalculator <InputImageType> ImageCalculatorType;

        if (m_TransformInitializationType != Identity)
        {
            InputImagePointer offsetRefImage = m_ReferenceImage;
            if (m_ReferenceMinimalValue < 0.0)
            {
                typedef itk::AddImageFilter <InputImageType, InputImageType, InputImageType> AddFilterType;
                typename AddFilterType::Pointer addFilter = AddFilterType::New();
                addFilter->SetInput(m_ReferenceImage);
                addFilter->SetConstant(- m_ReferenceMinimalValue);
                addFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
                addFilter->Update();

                offsetRefImage = addFilter->GetOutput();
                offsetRefImage->DisconnectPipeline();
            }

            typename ImageCalculatorType::Pointer fixedCalculator = ImageCalculatorType::New();
            fixedCalculator->SetImage(offsetRefImage);
            fixedCalculator->Compute();
            typename ImageCalculatorType::VectorType fixedBar = fixedCalculator->GetCenterOfGravity();

            InputImagePointer offsetFloatingImage = m_FloatingImage;
            if (m_ReferenceMinimalValue < 0.0)
            {
                typedef itk::AddImageFilter <InputImageType, InputImageType, InputImageType> AddFilterType;
                typename AddFilterType::Pointer addFilter = AddFilterType::New();
                addFilter->SetInput(m_FloatingImage);
                addFilter->SetConstant(- m_FloatingMinimalValue);
                addFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
                addFilter->Update();

                offsetFloatingImage = addFilter->GetOutput();
                offsetFloatingImage->DisconnectPipeline();
            }

            typename ImageCalculatorType::Pointer movingCalculator = ImageCalculatorType::New();
            movingCalculator->SetImage(offsetFloatingImage);
            movingCalculator->Compute();
            typename ImageCalculatorType::VectorType movingBar = movingCalculator->GetCenterOfGravity();

            m_InitialTransform->SetOffset(movingBar - fixedBar);

            if ((m_TransformInitializationType == ClosestTransform)&&(m_OutputTransformType != outTranslation))
            {
                typename ImageCalculatorType::VectorType fixedPrincipalMom = fixedCalculator->GetPrincipalMoments();
                typename ImageCalculatorType::ScalarType fixedMass = fixedCalculator->GetTotalMass();

                typename ImageCalculatorType::VectorType movingPrincipalMom = movingCalculator->GetPrincipalMoments();
                typename ImageCalculatorType::ScalarType movingMass = movingCalculator->GetTotalMass();

                vnl_matrix<double> scalMatrix(ImageDimension, ImageDimension, 0);
                itk::Vector<double, ImageDimension> scalOffset;
                double fixedScalFactor = 0;
                double movingScalFactor = 0;

                for (unsigned int i = 0; i < ImageDimension; ++i)
                {
                    fixedScalFactor += pow(fixedPrincipalMom[i] / fixedMass, 0.5);
                    movingScalFactor += pow(movingPrincipalMom[i] / movingMass, 0.5);
                }

                double scalingFactor = 1.0;
                if ((m_OutputTransformType == outAnisotropic_Sim)||(m_OutputTransformType == outAffine))
                    scalingFactor = movingScalFactor / fixedScalFactor;

                scalMatrix.fill_diagonal(scalingFactor);

                for (unsigned int i = 0; i < ImageDimension; ++i)
                    scalOffset[i] = movingBar[i] - scalingFactor * fixedBar[i];

                std::cout << scalMatrix << " " << scalOffset << std::endl;

                m_InitialTransform->SetMatrix(scalMatrix);
                m_InitialTransform->SetOffset(scalOffset);
            }
        }

        typename ResampleFilterType::Pointer tmpResample = ResampleFilterType::New();
        tmpResample->SetTransform(m_InitialTransform);
        tmpResample->SetInput(m_FloatingImage);

        tmpResample->SetSize(m_ReferenceImage->GetLargestPossibleRegion().GetSize());
        tmpResample->SetOutputOrigin(m_ReferenceImage->GetOrigin());
        tmpResample->SetOutputSpacing(m_ReferenceImage->GetSpacing());
        tmpResample->SetOutputDirection(m_ReferenceImage->GetDirection());
        tmpResample->SetDefaultPixelValue(m_FloatingMinimalValue);
        tmpResample->Update();

        initialFloatingImage = tmpResample->GetOutput();
        initialFloatingImage->DisconnectPipeline();

        using MinMaxImageFilterType = itk::MinimumMaximumImageFilter <InputImageType>;
        typename MinMaxImageFilterType::Pointer minMaxFilter = MinMaxImageFilterType::New();
        minMaxFilter->SetInput(initialFloatingImage);
        minMaxFilter->SetNumberOfWorkUnits(GetNumberOfWorkUnits());
        minMaxFilter->Update();

        if ((minMaxFilter->GetMinimum() == minMaxFilter->GetMaximum())&&(m_TransformInitializationType == Identity))
        {
            std::cout << "Identity initialization outputs an empty image, initializing with centers of mass" << std::endl;
            m_TransformInitializationType = GravityCenters;

            typename ImageCalculatorType::Pointer fixedCalculator = ImageCalculatorType::New();
            fixedCalculator->SetImage(m_ReferenceImage);
            fixedCalculator->Compute();
            typename ImageCalculatorType::VectorType fixedBar = fixedCalculator->GetCenterOfGravity();

            typename ImageCalculatorType::Pointer movingCalculator = ImageCalculatorType::New();
            movingCalculator->SetImage(m_FloatingImage);
            movingCalculator->Compute();
            typename ImageCalculatorType::VectorType movingBar = movingCalculator->GetCenterOfGravity();

            m_InitialTransform->SetOffset(movingBar - fixedBar);

            typename ResampleFilterType::Pointer tmpResample = ResampleFilterType::New();
            tmpResample->SetTransform(m_InitialTransform);
            tmpResample->SetInput(m_FloatingImage);

            tmpResample->SetSize(m_ReferenceImage->GetLargestPossibleRegion().GetSize());
            tmpResample->SetOutputOrigin(m_ReferenceImage->GetOrigin());
            tmpResample->SetOutputSpacing(m_ReferenceImage->GetSpacing());
            tmpResample->SetOutputDirection(m_ReferenceImage->GetDirection());
            tmpResample->SetDefaultPixelValue(m_FloatingMinimalValue);
            tmpResample->Update();

            initialFloatingImage = tmpResample->GetOutput();
            initialFloatingImage->DisconnectPipeline();
        }
    }

    // Create pyramid for Floating image
    m_ReferencePyramid = PyramidType::New();
    typename ResampleFilterType::Pointer refResampler = ResampleFilterType::New();

    bool invertInputs = (m_RegistrationPointLocation < 0.5) && (m_SymmetryType == Kissing);

    if (!invertInputs)
    {
        m_ReferencePyramid->SetInput(m_ReferenceImage);
        refResampler->SetDefaultPixelValue(m_ReferenceMinimalValue);
    }
    else
    {
        m_ReferencePyramid->SetInput(initialFloatingImage);
        refResampler->SetDefaultPixelValue(m_FloatingMinimalValue);
    }

    m_ReferencePyramid->SetNumberOfLevels(GetNumberOfPyramidLevels());
    m_ReferencePyramid->SetNumberOfWorkUnits(GetNumberOfWorkUnits());

    m_ReferencePyramid->SetImageResampler(refResampler);
    m_ReferencePyramid->Update();

    // Create pyramid for Floating image
    m_FloatingPyramid = PyramidType::New();
    typename ResampleFilterType::Pointer floResampler = ResampleFilterType::New();

    if (!invertInputs)
    {
        m_FloatingPyramid->SetInput(initialFloatingImage);
        floResampler->SetDefaultPixelValue(m_FloatingMinimalValue);
    }
    else
    {
        m_FloatingPyramid->SetInput(m_ReferenceImage);
        floResampler->SetDefaultPixelValue(m_ReferenceMinimalValue);
    }

    m_FloatingPyramid->SetNumberOfLevels(GetNumberOfPyramidLevels());
    m_FloatingPyramid->SetNumberOfWorkUnits(GetNumberOfWorkUnits());
    m_FloatingPyramid->SetImageResampler(floResampler);

    m_FloatingPyramid->Update();

    m_BlockGenerationPyramid = 0;
    if (m_BlockGenerationMask)
    {
        typedef anima::ResampleImageFilter<MaskImageType, MaskImageType,
                typename AgregatorType::ScalarType> MaskResampleFilterType;

        typename MaskResampleFilterType::Pointer maskResampler = MaskResampleFilterType::New();

        m_BlockGenerationPyramid = MaskPyramidType::New();
        m_BlockGenerationPyramid->SetImageResampler(maskResampler);
        m_BlockGenerationPyramid->SetInput(m_BlockGenerationMask);
        m_BlockGenerationPyramid->SetNumberOfLevels(GetNumberOfPyramidLevels());
        m_BlockGenerationPyramid->SetNumberOfWorkUnits(GetNumberOfWorkUnits());
        m_BlockGenerationPyramid->Update();
    }
}

} // end of namespace anima
