#pragma once
#include "animaBlockMatchingBaseRegistrationMethod.h"

#include <animaSmoothingRecursiveYvvGaussianImageFilter.h>

/* Scalar interpolators */
#include <animaResampleImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <animaFasterLinearInterpolateImageFunction.h>

/* Optimizers */
#include <animaVoxelExhaustiveOptimizer.h>
#include <animaBobyqaOptimizer.h>

/* Similarity measures */
#include <animaFastCorrelationImageToImageMetric.h>
#include <itkMeanSquaresImageToImageMetric.h>

/* Transforms */
#include <itkTranslationTransform.h>
#include <animaLogRigid3DTransform.h>
#include <animaSplitAffine3DTransform.h>

#include <animaMatrixLogExp.h>
#include <animaVelocityUtils.h>
#include <itkTimeProbe.h>

namespace anima
{

/**
 * Constructor
 */
template <typename TInputImage>
BlockMatchingBaseRegistrationMethod<TInputImage>
::BlockMatchingBaseRegistrationMethod()
{
    m_Abort = false;

    this->SetNumberOfRequiredOutputs( 1 );  // for the Transform

    m_FixedImage = 0; // has to be provided by the user.
    m_MovingImage = 0; // has to be provided by the user.

    m_InitialTransformParameters = ParametersType(1);
    m_InitialTransformParameters.Fill( 0.0f );

    m_AngleMax = M_PI;
    m_TranslateMax = 10;
    m_SkewMax = M_PI / 4.0;
    m_ScaleMax = 3;

    m_SearchRadius = 2;
    m_SearchAngleRadius = 5;
    m_SearchSkewRadius = 5;
    m_SearchScaleRadius = 0.1;

    m_MinimalTransformError = 0.0001;
    m_BlockScalarVarianceThreshold = 25;

    m_SVFElasticRegSigma = 0;

    TransformOutputPointer transformDecorator = static_cast <TransformOutputType *> (this->MakeOutput(0).GetPointer());
    this->itk::ProcessObject::SetNthOutput(0, transformDecorator.GetPointer());

    this->SetNumberOfThreads(this->GetMultiThreader()->GetNumberOfThreads());
}

template <typename TInputImage>
itk::ModifiedTimeType
BlockMatchingBaseRegistrationMethod<TInputImage>
::GetMTime() const
{
    unsigned long mtime = Superclass::GetMTime();
    unsigned long m;

    // Some of the following should be removed once ivars are put in the
    // input and output lists

    if (m_FixedImage)
    {
        m = m_FixedImage->GetMTime();
        mtime = (m > mtime ? m : mtime);
    }

    if (m_MovingImage)
    {
        m = m_MovingImage->GetMTime();
        mtime = (m > mtime ? m : mtime);
    }

    return mtime;
}

template <typename TInputImage>
typename BlockMatchingBaseRegistrationMethod<TInputImage>::BaseInputTransformType::Pointer
BlockMatchingBaseRegistrationMethod<TInputImage>
::GetBaseTransformPointer(unsigned int i)
{
    if (i < m_BaseTransformsPointers.size())
        return m_BaseTransformsPointers[i];

    return NULL;
}

template <typename TInputImage>
void
BlockMatchingBaseRegistrationMethod<TInputImage>
::SetBaseTransform(unsigned int i, BaseInputTransformType *trsf)
{
    if (i < m_BaseTransformsPointers.size())
        m_BaseTransformsPointers[i] = trsf;
}

/**
 * Initialize by setting the interconnects between components.
 */
template <typename TInputImage>
void
BlockMatchingBaseRegistrationMethod<TInputImage>
::Initialize() throw (itk::ExceptionObject)
{

    if( !m_FixedImage )
    {
        itkExceptionMacro(<<"FixedImage is not present");
    }

    if( !m_MovingImage )
    {
        itkExceptionMacro(<<"MovingImage is not present");
    }

    this->GetMultiThreader()->SetNumberOfThreads(this->GetNumberOfThreads());
}

template <typename TInputImage>
void
BlockMatchingBaseRegistrationMethod<TInputImage>
::SetupTransform(TransformPointer &optimizedTransform)
{
    if (m_Agregator->GetOutputTransformType() != AgregatorType::SVF)
    {
        if (m_InitialTransform)
        {
            optimizedTransform = AffineTransformType::New();
            optimizedTransform->SetParameters(m_InitialTransform->GetParameters());
        }
        else
        {
            typename AffineTransformType::Pointer tmpTrsf = AffineTransformType::New();
            tmpTrsf->SetIdentity();
            optimizedTransform = tmpTrsf;
        }
    }
    else
    {
        if (m_InitialTransform)
            optimizedTransform = m_InitialTransform;
        else
        {
            SVFTransformPointer tmpTrsf = SVFTransformType::New();
            tmpTrsf->SetIdentity();
            optimizedTransform = tmpTrsf;
        }
    }
}

template <typename TInputImage>
void
BlockMatchingBaseRegistrationMethod<TInputImage>
::SetBlockRegions(std::vector <ImageRegionType>& origins)
{
    m_BlockRegions = origins;

    unsigned int regionsNum = origins.size();

    m_BlockWeights.resize(regionsNum);
    m_BaseTransformsPointers.resize(regionsNum);
}

/**
 * Starts the Registration Process
 */
template <typename TInputImage>
void
BlockMatchingBaseRegistrationMethod<TInputImage>
::StartOptimization( void )
{
    switch (m_TransformKind)
    {
        case Translation:
            m_Agregator->SetInputTransformType(AgregatorType::TRANSLATION);
            break;
        case Rigid:
            m_Agregator->SetInputTransformType(AgregatorType::RIGID);
            break;
        case Affine:
            m_Agregator->SetInputTransformType(AgregatorType::AFFINE);
            break;
        default:
            std::cerr << "Input transform type not supported... Exiting..." << std::endl;
            break;
    }

    TransformPointer computedTransform = NULL;
    this->SetupTransform(computedTransform);

    this->GlobalParametersSetup();

    //progress management
    itk::ProgressReporter progress(this, 0, m_MaximumIterations);

    // Real work goes here
    InputImagePointer fixedResampled, movingResampled;
    for (unsigned int iterations = 0; iterations < m_MaximumIterations && !m_Abort; ++iterations)
    {
        // Resample fixed and moving image here
        this->ResampleImages(computedTransform, fixedResampled, movingResampled);

        // Perform one iteration of registration between the images
        // Calls pure virtual method that can use the block matching methods available here
        TransformPointer addOn;
        this->PerformOneIteration(fixedResampled, movingResampled, addOn);

        bool continueLoop = this->ComposeAddOnWithTransform(computedTransform,addOn);

        std::cout << "Iteration " << iterations << " done..." << std::endl;

        if (iterations != m_MaximumIterations - 1)
            progress.CompletedPixel();

        if (!continueLoop)
        {
            TransformOutputPointer transformDecorator = TransformOutputType::New();
            transformDecorator->Set(computedTransform.GetPointer());

            this->itk::ProcessObject::SetNthOutput(0, transformDecorator.GetPointer());

            return;
        }
    }

    if(m_Abort)
        m_StopConditionDescription << "Process aborted";
    else
        m_StopConditionDescription << "Maximum Iterations";

    TransformOutputPointer transformDecorator = TransformOutputType::New();
    transformDecorator->Set(computedTransform.GetPointer());

    this->itk::ProcessObject::SetNthOutput(0, transformDecorator.GetPointer());
}

template <typename TInputImage>
void
BlockMatchingBaseRegistrationMethod<TInputImage>
::ResampleImages(TransformType *currentTransform, InputImagePointer &refImage, InputImagePointer &movingImage)
{
    typedef anima::ResampleImageFilter <InputImageType,InputImageType,typename AgregatorType::ScalarType> ResampleFilterType;

    typename ResampleFilterType::Pointer resample = ResampleFilterType::New();

    resample->SetNumberOfThreads(this->GetNumberOfThreads());

    if (m_Agregator->GetOutputTransformType() == AgregatorType::SVF)
    {
        // Compute temporary field and set it to resampler
        DisplacementFieldTransformPointer dispTrsf = DisplacementFieldTransformType::New();
        SVFTransformType *svfCast = dynamic_cast<SVFTransformType *> (currentTransform);

        GetSVFExponential(svfCast,dispTrsf.GetPointer(),false);

        resample->SetTransform(dispTrsf);
    }
    else
        resample->SetTransform(currentTransform);

    resample->SetInput(m_MovingImage);

    resample->SetSize(m_FixedImage->GetLargestPossibleRegion().GetSize());
    resample->SetOutputOrigin(m_FixedImage->GetOrigin());
    resample->SetOutputSpacing(m_FixedImage->GetSpacing());
    resample->SetOutputDirection(m_FixedImage->GetDirection());
    resample->SetDefaultPixelValue(0);

    resample->Update();

    movingImage = resample->GetOutput();
    movingImage->DisconnectPipeline();

    typename ResampleFilterType::Pointer resampleFixed = ResampleFilterType::New();

    resampleFixed->SetNumberOfThreads(this->GetNumberOfThreads());

    // Compute temporary field and set it to resampler
    if (m_Agregator->GetOutputTransformType() == AgregatorType::SVF)
    {
        DisplacementFieldTransformPointer dispTrsf = DisplacementFieldTransformType::New();
        SVFTransformType *svfCast = dynamic_cast<SVFTransformType *> (currentTransform);

        GetSVFExponential(svfCast,dispTrsf.GetPointer(),true);

        resampleFixed->SetTransform(dispTrsf);
    }
    else
    {
        AffineTransformType *affCast = dynamic_cast<AffineTransformType *> (currentTransform);
        TransformPointer reverseTrsf = affCast->GetInverseTransform();
        resampleFixed->SetTransform(reverseTrsf);
    }

    resampleFixed->SetInput(m_FixedImage);

    resampleFixed->SetSize(m_FixedImage->GetLargestPossibleRegion().GetSize());
    resampleFixed->SetOutputOrigin(m_FixedImage->GetOrigin());
    resampleFixed->SetOutputSpacing(m_FixedImage->GetSpacing());
    resampleFixed->SetOutputDirection(m_FixedImage->GetDirection());
    resampleFixed->SetDefaultPixelValue(0);

    resampleFixed->Update();

    refImage = resampleFixed->GetOutput();
    refImage->DisconnectPipeline();
}

template <typename TInputImage>
bool
BlockMatchingBaseRegistrationMethod<TInputImage>
::ComposeAddOnWithTransform(TransformType *computedTransform, TransformType *addOn)
{
    if (m_Agregator->GetOutputTransformType() != AgregatorType::SVF)
    {
        typename TransformType::ParametersType oldPars = computedTransform->GetParameters();

        AffineTransformType *tmpTrsf = dynamic_cast<AffineTransformType *>(computedTransform);
        AffineTransformType *tmpAddOn = dynamic_cast<AffineTransformType *>(addOn);
        tmpTrsf->Compose(tmpAddOn, true);

        typename TransformType::ParametersType newPars = computedTransform->GetParameters();

        // Compute the distance between consecutive solutions, until a certain threshold
        double err = 0;
        for (unsigned int i = 0; i < newPars.Size(); ++i)
            err += pow(newPars[i] - oldPars[i], 2.);

        if (err < m_MinimalTransformError)
        {
            m_StopConditionDescription << "Minimal Transform Error Reached " << err;
            return false;
        }
    }
    else
    {
        // Add update to current velocity field (cf. Vercauteren et al, 2008)
        SVFTransformType *tmpTrsf = dynamic_cast<SVFTransformType *>(computedTransform);
        SVFTransformType *tmpAddOn = dynamic_cast<SVFTransformType *>(addOn);

        anima::composeSVF(tmpTrsf,tmpAddOn);
        if (m_SVFElasticRegSigma > 0)
        {
            typedef typename SVFTransformType::VectorFieldType VelocityFieldType;

            typedef anima::SmoothingRecursiveYvvGaussianImageFilter <VelocityFieldType, VelocityFieldType> SmoothingFilterType;
            typename SmoothingFilterType::Pointer smootherPtr = SmoothingFilterType::New();

            smootherPtr->SetInput(tmpTrsf->GetParametersAsVectorField());
            smootherPtr->SetSigma(m_SVFElasticRegSigma);
            smootherPtr->SetNumberOfThreads(this->GetNumberOfThreads());

            smootherPtr->Update();

            typename VelocityFieldType::Pointer tmpSmoothed = smootherPtr->GetOutput();
            tmpSmoothed->DisconnectPipeline();

            tmpTrsf->SetParametersAsVectorField(tmpSmoothed);
        }
    }

    return true;
}

/**
 * Threaded block match
 */
template <typename TInputImage>
ITK_THREAD_RETURN_TYPE
BlockMatchingBaseRegistrationMethod<TInputImage>
::ThreadedMatching(void *arg)
{
    itk::MultiThreader::ThreadInfoStruct *threadArgs = (itk::MultiThreader::ThreadInfoStruct *)arg;
    // Set the registration procedure
    ThreadedMatchData* data = (ThreadedMatchData *)threadArgs->UserData;

    data->BlockMatch->BlockMatch(threadArgs->ThreadID, threadArgs->NumberOfThreads, data->fixedImage, data->movingImage);

    return NULL;
}

template <typename TInputImage>
void
BlockMatchingBaseRegistrationMethod<TInputImage>
::GlobalParametersSetup()
{
    switch (m_TransformKind)
    {
        case Translation:
        {
            typedef itk::TranslationTransform <double, InputImageType::ImageDimension> itkTransformType;
            typename itkTransformType::Pointer tr = itkTransformType::New();
            tr->SetIdentity();

            m_InitialTransformParameters = tr->GetParameters();

            break;
        }

        case Rigid:
        {
            if (InputImageType::ImageDimension == 3)
            {
                typedef anima::LogRigid3DTransform<double> itkTransformType;
                typename itkTransformType::Pointer tr = itkTransformType::New();
                tr->SetIdentity();

                m_InitialTransformParameters = tr->GetParameters();
            }
            else
                std::cout << "Only Rigid 3D transforms handled right now." << std::endl;
            break;
        }

        case Affine:
        {
            typedef anima::SplitAffine3DTransform <double> itkTransformType;
            typename itkTransformType::Pointer tr = itkTransformType::New();
            tr->SetIdentity();

            m_InitialTransformParameters = tr->GetParameters();

            break;
        }

    }

    switch (m_MetricKind)
    {
        case MeanSquares:
            m_MaximizeMetric = false;
            break;

        case Correlation:
        case SquaredCorrelation:
            m_MaximizeMetric = true;
            break;
    }
}

template <typename TInputImage>
void
BlockMatchingBaseRegistrationMethod<TInputImage>
::SetupBlockParameters(InterpolatorPointer& interpolator, MetricPointer &metric, OptimizerPointer &optimizer,
                       unsigned startBlock, unsigned endBlock)

{
    // Transformation types handled, kind of redundant with itk transforms but much easier to use after
    // TRANSLATION -> itk::TranslationTransform
    // RIGID -> itk::LogRigid3DTransform
    // AFFINE -> itk::SplitAffine3DTransform

    switch (m_TransformKind)
    {
        case Translation:
        {
            typedef itk::TranslationTransform <double, InputImageType::ImageDimension> itkTransformType;
            typename itkTransformType::Pointer tr = itkTransformType::New();

            tr->SetIdentity();

            for (unsigned i = startBlock; i < endBlock; i++)
            {
                m_BaseTransformsPointers[i] = itkTransformType::New();
                m_BaseTransformsPointers[i]->SetParameters(tr->GetParameters());
            }

            break;
        }

        case Rigid:
        {
            if (InputImageType::ImageDimension == 3)
            {
                typedef anima::LogRigid3DTransform<double> itkTransformType;
                typename itkTransformType::Pointer tr = itkTransformType::New();

                tr->SetIdentity();
                itk::ContinuousIndex <double,InputImageType::ImageDimension> tmpIndex;

                for (unsigned i = startBlock; i < endBlock; i++)
                {
                    typename itkTransformType::Pointer tmpTr = itkTransformType::New();
                    tmpTr->SetCenter(this->m_Agregator->GetInputOrigin(i));

                    m_BaseTransformsPointers[i] = tmpTr;
                    m_BaseTransformsPointers[i]->SetParameters(tr->GetParameters());
                }

            }
            else
                std::cout << "Only Rigid 3D transforms handled right now." << std::endl;
            break;
        }

        case Affine:
        {
            if (InputImageType::ImageDimension == 3)
            {
                typedef anima::SplitAffine3DTransform <double> itkTransformType;
                typename itkTransformType::Pointer tr = itkTransformType::New();

                tr->SetIdentity();
                itk::ContinuousIndex <double,InputImageType::ImageDimension> tmpIndex;

                for (unsigned i = startBlock; i < endBlock; i++)
                {
                    typename itkTransformType::Pointer tmpTr = itkTransformType::New();
                    tmpTr->SetCenter(m_Agregator->GetInputOrigin(i));

                    m_BaseTransformsPointers[i] = tmpTr;
                    m_BaseTransformsPointers[i]->SetParameters(tr->GetParameters());
                }
            }
            else
                std::cout << "Only 3D affine transforms handled right now." << std::endl;

            break;
        }

    }

    switch (m_MetricKind)
    {
        case Correlation:
        case SquaredCorrelation:
        {
            typedef anima::FastCorrelationImageToImageMetric <InputImageType,InputImageType> MetricType;

            typename MetricType::Pointer tmpMetric = MetricType::New();
            tmpMetric->SetSquaredCorrelation(m_MetricKind == SquaredCorrelation);
            tmpMetric->SetVarianceThreshold(this->GetBlockScalarVarianceThreshold());

            metric = tmpMetric;
            break;

        }

        case MeanSquares:
        {
            typedef itk::MeanSquaresImageToImageMetric <InputImageType,InputImageType> MetricType;

            metric = MetricType::New();
            break;
        }

        default:
            itkExceptionMacro("Metric not implemented yet");
            break;
    }

    typedef anima::FasterLinearInterpolateImageFunction<InputImageType,double> LocalInterpolatorType;
    interpolator = LocalInterpolatorType::New();

    metric->SetInterpolator(interpolator);

    switch (m_OptimizerKind)
    {
        case Bobyqa:
        {
            typedef anima::BobyqaOptimizer LocalOptimizerType;
            optimizer = LocalOptimizerType::New();
            LocalOptimizerType * tmpOpt = (LocalOptimizerType *)optimizer.GetPointer();
            tmpOpt->SetRhoBegin(m_SearchRadius);
            tmpOpt->SetRhoEnd(m_FinalRadius);

            tmpOpt->SetNumberSamplingPoints(m_BaseTransformsPointers[startBlock]->GetNumberOfParameters() + 2);
            // Maximum sampling of cost function ((Par+1)(Par+2))/2
            tmpOpt->SetMaximumIteration(m_OptimizerMaximumIterations);
            typename InputImageType::SpacingType fixedSpacing = m_FixedImage->GetSpacing();

            LocalOptimizerType::ScalesType tmpScales(m_BaseTransformsPointers[startBlock]->GetNumberOfParameters());
            LocalOptimizerType::ScalesType lowerBounds(m_BaseTransformsPointers[startBlock]->GetNumberOfParameters());
            LocalOptimizerType::ScalesType upperBounds(m_BaseTransformsPointers[startBlock]->GetNumberOfParameters());

            if (m_TransformKind == Rigid)
            {
                for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
                {
                    tmpScales[i] = m_SearchRadius * 180.0 / (m_SearchAngleRadius * M_PI);
                    lowerBounds[i] = - m_AngleMax * M_PI / 180.0;
                    upperBounds[i] = m_AngleMax * M_PI / 180.0;

                    tmpScales[i+InputImageType::ImageDimension] = 1.0 / fixedSpacing[i];

                    lowerBounds[i+InputImageType::ImageDimension] = - m_TranslateMax * fixedSpacing[i];
                    upperBounds[i+InputImageType::ImageDimension] = m_TranslateMax * fixedSpacing[i];
                }
            }
            else if (m_TransformKind == Translation)
            {
                for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
                {
                    tmpScales[i] = 1.0 / fixedSpacing[i];

                    lowerBounds[i] = - m_TranslateMax * fixedSpacing[i];
                    upperBounds[i] = m_TranslateMax * fixedSpacing[i];
                }
            }
            else // Affine
            {
                // There are 12 parameters: 3 angles, 3 translations, 3 scales, 3 skew factors
                for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
                {
                    // Angles
                    tmpScales[i] = m_SearchRadius * 180.0 / (m_SearchAngleRadius * M_PI);
                    lowerBounds[i] = - m_AngleMax * M_PI / 180.0;
                    upperBounds[i] = m_AngleMax * M_PI / 180.0;

                    // Translations
                    tmpScales[InputImageType::ImageDimension + i] = 1.0 / fixedSpacing[i];

                    lowerBounds[InputImageType::ImageDimension + i] = - m_TranslateMax * fixedSpacing[i];
                    upperBounds[InputImageType::ImageDimension + i] = m_TranslateMax * fixedSpacing[i];

                    // Scales
                    tmpScales[2 * InputImageType::ImageDimension + i] = m_SearchRadius / m_SearchScaleRadius;

                    lowerBounds[2 * InputImageType::ImageDimension + i] = - log (m_ScaleMax);
                    upperBounds[2 * InputImageType::ImageDimension + i] = log (m_ScaleMax);

                    // Skews
                    tmpScales[3 * InputImageType::ImageDimension + i] = m_SearchRadius * 180.0 / (m_SearchSkewRadius * M_PI);

                    lowerBounds[3 * InputImageType::ImageDimension + i] = - m_SkewMax * M_PI / 180.0;
                    upperBounds[3 * InputImageType::ImageDimension + i] = m_SkewMax * M_PI / 180.0;
                }
            }

            tmpOpt->SetScales(tmpScales);
            tmpOpt->SetLowerBounds(lowerBounds);
            tmpOpt->SetUpperBounds(upperBounds);
            tmpOpt->SetMaximize(m_MaximizeMetric);
            metric->ComputeGradientOff();
            break;
        }

        case Exhaustive:
        {
            typedef anima::VoxelExhaustiveOptimizer LocalOptimizerType;
            optimizer = LocalOptimizerType::New();

            LocalOptimizerType * tmpOpt = (LocalOptimizerType *)optimizer.GetPointer();
            LocalOptimizerType::StepsType steps(m_BaseTransformsPointers[startBlock]->GetNumberOfParameters());

            typename InputImageType::SpacingType fixedSpacing = m_FixedImage->GetSpacing();
            typename InputImageType::DirectionType fixedDirection = m_FixedImage->GetDirection();

            vnl_matrix <double> geometry(InputImageType::ImageDimension,InputImageType::ImageDimension);
            for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
            {
                double tmpVoxSize = fixedSpacing[i];

                for (unsigned int j = 0;j < InputImageType::ImageDimension;++j)
                    geometry(j,i) = tmpVoxSize*fixedDirection(j,i);
            }

            tmpOpt->SetGeometry(geometry);

            LocalOptimizerType::ScalesType tmpScales(InputImageType::ImageDimension);

            for (unsigned i = 0; i < steps.Size(); ++i)
            {
                steps[i] = m_SearchRadius;
                tmpScales[i] = m_StepSize;
            }

            tmpOpt->SetGeometry(m_Geometry);
            tmpOpt->SetNumberOfSteps(steps);
            tmpOpt->SetScales (tmpScales);
            tmpOpt->SetMaximize(m_MaximizeMetric);

            metric->ComputeGradientOff();
            break;
        }
    }
}

template <typename TInputImage>
void
BlockMatchingBaseRegistrationMethod<TInputImage>
::AdditionalBlockMatchingSetup(MetricPointer &metric, unsigned int blockNum)
{
    if (m_MetricKind != MeanSquares)
        ((anima::FastCorrelationImageToImageMetric<InputImageType, InputImageType> *)metric.GetPointer())->PreComputeFixedValues();
}

template <typename TInputImage>
void
BlockMatchingBaseRegistrationMethod<TInputImage>
::BlockMatch(unsigned threadId, unsigned NbThreads, InputImageType *fixedImage,
             InputImageType *movingImage)
{
    InterpolatorPointer interpolator;
    MetricPointer metric;
    OptimizerPointer optimizer;

    unsigned totalRegionSize = m_BlockRegions.size();
    unsigned step = totalRegionSize / NbThreads;

    unsigned startBlock = threadId*step;
    unsigned endBlock = (1+threadId)*step;

    if (threadId+1==NbThreads) // Ensure we got up to the last block
        endBlock = totalRegionSize;

    this->SetupBlockParameters(interpolator, metric, optimizer, startBlock, endBlock);

    metric->SetFixedImage(fixedImage);
    metric->SetMovingImage(movingImage);
    interpolator->SetInputImage(movingImage);

    itkDebugMacro(<<"Threaded Block Match " << threadId << " / " << NbThreads << " starting @ " << startBlock << " up to " << endBlock );

    // Loop over the desired blocks
    for (unsigned int block = startBlock; block < endBlock; ++block)
    {

        metric->SetFixedImageRegion(m_BlockRegions[block]);
        metric->SetTransform(m_BaseTransformsPointers[block]);
        metric->Initialize();

        this->AdditionalBlockMatchingSetup(metric, block);

        optimizer->SetCostFunction(metric);
        optimizer->SetInitialPosition(m_InitialTransformParameters);

        try
        {
            optimizer->StartOptimization();
        }
        catch( itk::ExceptionObject & err )
        {
            std::cout << "ExceptionObject caught Image Block Match! " << err << std::endl;
            m_BlockWeights[block] = 0;
            continue;
        }

        m_BaseTransformsPointers[block]->SetParameters(optimizer->GetCurrentPosition());

        if (m_MetricKind != MeanSquares)
        {
            double val = 0;
            switch (m_OptimizerKind)
            {
                case Bobyqa:
                    val = ((anima::BobyqaOptimizer *)optimizer.GetPointer())->GetValue();
                    break;

                case Exhaustive:
                default:
                    if (m_MaximizeMetric)
                        val = ((anima::VoxelExhaustiveOptimizer *)optimizer.GetPointer())->GetCurrentCost();

                    break;
            }

            m_BlockWeights[block] =  (val >= 0) ? val : -val;
        }
        else
            m_BlockWeights[block] =  1;

        itkDebugMacro(<<  m_BaseTransformsPointers[block]->GetParameters() << "  " << m_BlockWeights[block] << std::endl);
    }
}

/**
 * PrintSelf
 */
template <typename TInputImage>
void
BlockMatchingBaseRegistrationMethod<TInputImage>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
    Superclass::PrintSelf( os, indent );
    os << indent << "Optimizer: " << (char*)(this->m_OptimizerKind == Self::Bobyqa ? "Bobyqa" : "Exhaustive") << std::endl;
    os << indent << "Transform: " << (char*)(this->m_TransformKind == Self::Translation ? "Translations" : this->m_TransformKind == Self::Rigid ? "Rigid" : "Affine") << std::endl;

    os << indent << "Fixed Image: " << m_FixedImage.GetPointer() << std::endl;
    os << indent << "Moving Image: " << m_MovingImage.GetPointer() << std::endl;

    os << indent << "Search Radius: " << m_SearchRadius << std::endl;
    os << indent << "Final Radius: " << m_FinalRadius << std::endl;
    os << indent << "Step Size : " << m_StepSize << std::endl;

    os << indent << "Maximum Iterations: " << m_MaximumIterations << std::endl;

    os << indent << "Optimizer Maximum Iterations: " << m_OptimizerMaximumIterations << std::endl;
    os << indent << "Initial Transform Parameters: " << m_InitialTransformParameters << std::endl;
}

/*
 * Generate Data
 */
template <typename TInputImage>
void
BlockMatchingBaseRegistrationMethod<TInputImage>
::GenerateData()
{
    m_Abort = false;

    ParametersType empty(1);
    empty.Fill( 0.0 );
    try
    {
        // initialize the interconnects between components
        this->Initialize();

    }
    catch (itk::ExceptionObject &err)
    {
        // pass exception to caller
        throw err;
    }

    this->StartOptimization();
}

/**
 *  Get Output
 */
template <typename TInputImage>
typename BlockMatchingBaseRegistrationMethod<TInputImage>::TransformOutputType *
BlockMatchingBaseRegistrationMethod<TInputImage>
::GetOutput()
{
    return static_cast <TransformOutputType *> (this->ProcessObject::GetOutput(0));
}

template <typename TInputImage>
itk::DataObject::Pointer
BlockMatchingBaseRegistrationMethod<TInputImage>
::MakeOutput(DataObjectPointerArraySizeType output)
{
    switch (output)
    {
        case 0:
            return static_cast <itk::DataObject*> (TransformOutputType::New().GetPointer());
            break;
        default:
            itkExceptionMacro("MakeOutput request for an output number larger than the expected number of outputs");
            return 0;
    }
}

template <typename TInputImage>
void
BlockMatchingBaseRegistrationMethod<TInputImage>
::SetFixedImage (InputImageType *fixedImage)
{
    itkDebugMacro("Setting Fixed Image to " << fixedImage);

    if (m_FixedImage.GetPointer() != fixedImage)
    {
        m_FixedImage = fixedImage;
        m_FixedImage->DisconnectPipeline();

        // Reset geometry info
        m_Geometry.set_size(TInputImage::ImageDimension,TInputImage::ImageDimension);
        m_Geometry.set_identity();

        for (unsigned int j = 0;j < TInputImage::ImageDimension;++j)
        {
            double spacing = m_FixedImage->GetSpacing()[j];
            for (unsigned int i = 0;i < TInputImage::ImageDimension;++i)
                m_Geometry(j,i) = spacing * m_FixedImage->GetDirection()(j,i);
        }

        this->Modified();
    }
}

template <typename TInputImage>
void
BlockMatchingBaseRegistrationMethod<TInputImage>
::SetMovingImage (InputImageType * movingImage)
{
    itkDebugMacro("Setting Moving Image to " << movingImage);

    if (m_MovingImage.GetPointer() != movingImage)
    {
        m_MovingImage = movingImage;
        m_MovingImage->DisconnectPipeline();

        this->Modified();
    }
}

template <typename TInputImage>
const std::string
BlockMatchingBaseRegistrationMethod<TInputImage>
::GetStopConditionDescription() const
{
    return m_StopConditionDescription.str();
}

} // end namespace anima
