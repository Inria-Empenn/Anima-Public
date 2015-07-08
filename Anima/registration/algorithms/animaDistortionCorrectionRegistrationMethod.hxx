#pragma once
#include "animaDistortionCorrectionRegistrationMethod.h"

#include <itkImageFileWriter.h>

#include <animaSmoothingRecursiveYvvGaussianImageFilter.h>
#include <itkComposeDisplacementFieldsImageFilter.h>

#include <itkSubtractImageFilter.h>
#include <itkMultiplyImageFilter.h>

/* Similarity measures */
#include <animaFastCorrelationImageToImageMetric.h>
#include <animaFastMeanSquaresImageToImageMetric.h>

/* Scalar interpolators */
#include <animaResampleImageFilter.h>
#include <animaFasterLinearInterpolateImageFunction.h>

/* Optimizers */
#include <animaVoxelExhaustiveOptimizer.h>
#include <animaNewuoaOptimizer.h>
#include <animaBobyqaOptimizer.h>

/* Transforms */
#include <animaDirectionScaleSkewTransform.h>

/* Agregators */
#include <animaBalooSVFTransformAgregator.h>
#include <animaDenseSVFTransformAgregator.h>

#include <animaBlockMatchInitializer.h>
#include <animaMatrixLogExp.h>
#include <animaVelocityUtils.h>
#include <itkTimeProbe.h>

#include <itkCommand.h>
#include <itkImageFileWriter.h>
#include <vnl/vnl_trace.h>

class CommandIterationUpdate : public itk::Command
{
public:
    typedef  CommandIterationUpdate   Self;
    typedef  itk::Command             Superclass;
    typedef itk::SmartPointer<Self>   Pointer;
    itkNewMacro( Self );
protected:
    CommandIterationUpdate() {};
public:
    typedef   itk::SingleValuedNonLinearOptimizer         OptimizerType;
    typedef   const OptimizerType *                  OptimizerPointer;
    
    void Execute(itk::Object *caller, const itk::EventObject & event)
    {
        Execute( (const itk::Object *)caller, event);
    }
    
    void Execute(const itk::Object * object, const itk::EventObject & event)
    {
        OptimizerPointer optimizer =
                dynamic_cast< OptimizerPointer >( object );
        if( ! itk::IterationEvent().CheckEvent( &event ) )
        {
            return;
        }
        
        //      std::cout << optimizer->GetCurrentIteration() << "   ";
        std::cout << "\t\t"<< optimizer->GetValue(optimizer->GetCurrentPosition() ) << "   ";
        std::cout << optimizer->GetCurrentPosition() << std::endl;
    }
};

namespace anima
{
/**
 * Constructor
 */
template <typename TInputImage>
DistortionCorrectionRegistrationMethod<TInputImage>
::DistortionCorrectionRegistrationMethod()
{
    m_Abort = false;

    this->SetNumberOfRequiredOutputs( 1 );  // for the Transform

    m_BackwardImage   = 0; // has to be provided by the user.
    m_ForwardImage  = 0; // has to be provided by the user.

    m_InitialTransformParameters = ParametersType(1);
    m_LastTransformParameters = ParametersType(1);

    m_InitialTransformParameters.Fill( 0.0f );
    m_LastTransformParameters.Fill( 0.0f );

    m_TranslateMax = 10;
    m_ScaleMax = std::log(5);
    m_SkewMax = std::tan(M_PI / 3.0);

    m_SearchRadius = 2;
    m_SearchScaleRadius = 0.1;
    m_SearchSkewRadius = 0.1;

    m_SVFElasticRegSigma = 0;

    m_TransformationKind = DirectionScaleSkew;
    m_MetricKind = SquaredCorrelation;
    m_MaximiseMetric = true;

    m_BlockPercentageKept = 0.8;
    m_BlockSize = 5;
    m_BlockSpacing = 2;
    m_BlockScalarVarianceThreshold = 0;

    m_UseTransformationDam = true;
    m_DamDistance = 6;

    TransformOutputPointer transformDecorator = static_cast< TransformOutputType * >(this->MakeOutput(0).GetPointer());
    this->ProcessObject::SetNthOutput(0, transformDecorator.GetPointer());

    this->SetNumberOfThreads(this->GetMultiThreader()->GetNumberOfThreads());
}

template <typename TInputImage>
itk::ModifiedTimeType
DistortionCorrectionRegistrationMethod<TInputImage>
::GetMTime() const
{
    unsigned long mtime = Superclass::GetMTime();
    unsigned long m;

    // Some of the following should be removed once ivars are put in the
    // input and output lists

    if (m_BackwardImage)
    {
        m = m_BackwardImage->GetMTime();
        mtime = (m > mtime ? m : mtime);
    }

    if (m_ForwardImage)
    {
        m = m_ForwardImage->GetMTime();
        mtime = (m > mtime ? m : mtime);
    }

    return mtime;

}

/**
     * Initialize by setting the interconnects between components.
     */
template <typename TInputImage>
void
DistortionCorrectionRegistrationMethod<TInputImage>
::Initialize() throw (itk::ExceptionObject)
{

    if( !m_BackwardImage )
    {
        itkExceptionMacro(<<"BackwardImage is not present");
    }

    if( !m_ForwardImage )
    {
        itkExceptionMacro(<<"ForwardImage is not present");
    }

    // TO DO : test if data is in the same frame, if not complain
    this->GetMultiThreader()->SetNumberOfThreads(this->GetNumberOfThreads());
}

/**
     * Starts the Registration Process
     */
template <typename TInputImage>
void
DistortionCorrectionRegistrationMethod<TInputImage>
::StartOptimization( void )
{
    typedef typename DisplacementFieldTransformType::VectorFieldType VectorFieldType;

    std::fill(this->m_BlockWeights.begin(), this->m_BlockWeights.end(), 1.);

    // TODO : when doing direction, this should be translation...
    this->m_Agregator->SetInputTransformType(AgregatorType::AFFINE);

    DisplacementFieldTransformPointer computedTransform = NULL;
    if (m_CurrentTransform)
        computedTransform = m_CurrentTransform;
    else
    {
        computedTransform = DisplacementFieldTransformType::New();
        computedTransform->SetIdentity();
    }

    this->GlobalParametersSetup();

    //progress management
    itk::ProgressReporter progress(this, 0, this->m_MaximumIterations);

    // Real work goes here
    InputImagePointer backwardResampled, forwardResampled;
    for (unsigned int iterations = 0; iterations < this->m_MaximumIterations && !m_Abort; ++iterations)
    {
        // Resample fixed and moving image here
        this->ResampleImages(computedTransform, backwardResampled, forwardResampled);

        // Perform one registration step moving -> fixed
        SVFTransformPointer positiveAddOn = SVFTransformType::New();

        this->PerformOneIteration(forwardResampled, backwardResampled, positiveAddOn);

        // Now compute positive and negative updated transform
        DisplacementFieldTransformPointer positiveDispTrsf = DisplacementFieldTransformType::New();
        GetSVFExponential(positiveAddOn.GetPointer(),positiveDispTrsf.GetPointer(),false);

        DisplacementFieldTransformPointer negativeDispTrsf = DisplacementFieldTransformType::New();
        GetSVFExponential(positiveAddOn.GetPointer(),negativeDispTrsf.GetPointer(),true);

        anima::composeDistortionCorrections<typename AgregatorType::ScalarType, InputImageType::ImageDimension>
                (computedTransform,positiveDispTrsf,negativeDispTrsf,this->GetNumberOfThreads());

        // Smooth (elastic)
        if (m_SVFElasticRegSigma > 0)
        {
            typedef anima::SmoothingRecursiveYvvGaussianImageFilter <VectorFieldType, VectorFieldType> SmoothingFilterType;
            typename SmoothingFilterType::Pointer smootherPtr = SmoothingFilterType::New();

            smootherPtr->SetInput(computedTransform->GetParametersAsVectorField());
            smootherPtr->SetSigma(m_SVFElasticRegSigma);
            smootherPtr->SetNumberOfThreads(this->GetNumberOfThreads());

            smootherPtr->Update();

            typename VectorFieldType::Pointer tmpSmoothed = smootherPtr->GetOutput();
            tmpSmoothed->DisconnectPipeline();

            computedTransform->SetParametersAsVectorField(tmpSmoothed);
        }

        std::cout << "Iteration " << iterations << " done..." << std::endl;

        if (iterations != m_MaximumIterations - 1)
            progress.CompletedPixel();
    }

    if(m_Abort)
        m_StopConditionDescription << "Process aborted";
    else m_StopConditionDescription << "Maximum Iterations";

    this->GetOutput()->Set(computedTransform);
}

template <typename TInputImage>
void
DistortionCorrectionRegistrationMethod<TInputImage>
::ResampleImages(TransformType *currentTransform, InputImagePointer &refImage, InputImagePointer &movingImage)
{
    typedef typename DisplacementFieldTransformType::VectorFieldType VectorFieldType;
    typedef itk::MultiplyImageFilter <VectorFieldType,itk::Image <float, TInputImage::ImageDimension>, VectorFieldType> MultiplyFilterType;
    typedef itk::ComposeDisplacementFieldsImageFilter <VectorFieldType,VectorFieldType> ComposeFilterType;
    typedef typename itk::ImageRegionIterator <VectorFieldType> VectorFieldIterator;
    typedef typename VectorFieldType::PixelType VectorType;

    DisplacementFieldTransformPointer positiveTrsf;
    DisplacementFieldTransformPointer negativeTrsf;

    if (m_InitialTransform)
    {
        // Here compose and make things straight again
        positiveTrsf = DisplacementFieldTransformType::New();

        typename ComposeFilterType::Pointer composePositiveFilter = ComposeFilterType::New();
        composePositiveFilter->SetWarpingField(currentTransform->GetParametersAsVectorField());
        composePositiveFilter->SetDisplacementField(m_InitialTransform->GetParametersAsVectorField());
        composePositiveFilter->SetNumberOfThreads(this->GetNumberOfThreads());

        composePositiveFilter->Update();
        positiveTrsf->SetParametersAsVectorField(composePositiveFilter->GetOutput());

        typename MultiplyFilterType::Pointer multiplyInitFilter = MultiplyFilterType::New();
        multiplyInitFilter->SetInput(m_InitialTransform->GetParametersAsVectorField());
        multiplyInitFilter->SetConstant(-1.0);
        multiplyInitFilter->SetNumberOfThreads(this->GetNumberOfThreads());

        multiplyInitFilter->Update();

        typename MultiplyFilterType::Pointer multiplyCurrentFilter = MultiplyFilterType::New();
        multiplyCurrentFilter->SetInput(currentTransform->GetParametersAsVectorField());
        multiplyCurrentFilter->SetConstant(-1.0);
        multiplyCurrentFilter->SetNumberOfThreads(this->GetNumberOfThreads());

        multiplyCurrentFilter->Update();

        typename ComposeFilterType::Pointer composeNegativeFilter = ComposeFilterType::New();
        composeNegativeFilter->SetWarpingField(multiplyCurrentFilter->GetOutput());
        composeNegativeFilter->SetDisplacementField(multiplyInitFilter->GetOutput());
        composeNegativeFilter->SetNumberOfThreads(this->GetNumberOfThreads());

        composeNegativeFilter->Update();
        negativeTrsf = DisplacementFieldTransformType::New();
        negativeTrsf->SetParametersAsVectorField(composeNegativeFilter->GetOutput());

        VectorFieldIterator positiveItr(const_cast <VectorFieldType *> (positiveTrsf->GetParametersAsVectorField()),
                                        positiveTrsf->GetParametersAsVectorField()->GetLargestPossibleRegion());

        VectorFieldIterator negativeItr(const_cast <VectorFieldType *> (negativeTrsf->GetParametersAsVectorField()),
                                        negativeTrsf->GetParametersAsVectorField()->GetLargestPossibleRegion());

        // And compose them to get a transformation conform to disto correction requirements
        VectorType tmpVec;
        while (!positiveItr.IsAtEnd())
        {
            tmpVec = 0.5 * (positiveItr.Get() - negativeItr.Get());
            positiveItr.Set(tmpVec);
            negativeItr.Set(- tmpVec);

            ++positiveItr;
            ++negativeItr;
        }
    }
    else
    {
        // Just use the current transform which is already in good shape
        positiveTrsf = currentTransform;
        negativeTrsf = DisplacementFieldTransformType::New();

        typename MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
        multiplyFilter->SetInput(positiveTrsf->GetParametersAsVectorField());
        multiplyFilter->SetConstant(-1.0);
        multiplyFilter->SetNumberOfThreads(this->GetNumberOfThreads());

        multiplyFilter->Update();
        negativeTrsf->SetParametersAsVectorField(multiplyFilter->GetOutput());
    }


    typedef anima::ResampleImageFilter <InputImageType,InputImageType,
            typename AgregatorType::ScalarType> ResampleFilterType;

    typename ResampleFilterType::Pointer resample = ResampleFilterType::New();

    resample->SetNumberOfThreads(this->GetNumberOfThreads());

    resample->SetTransform(positiveTrsf);
    resample->SetInput(m_ForwardImage);

    resample->SetSize(m_BackwardImage->GetLargestPossibleRegion().GetSize());
    resample->SetOutputOrigin(m_BackwardImage->GetOrigin());
    resample->SetOutputSpacing(m_BackwardImage->GetSpacing());
    resample->SetOutputDirection(m_BackwardImage->GetDirection());
    resample->SetDefaultPixelValue(0);
    resample->SetScaleIntensitiesWithJacobian(true);

    resample->Update();

    movingImage = resample->GetOutput();
    movingImage->DisconnectPipeline();

    typename ResampleFilterType::Pointer resampleFixed = ResampleFilterType::New();

    resampleFixed->SetNumberOfThreads(this->GetNumberOfThreads());

    // Compute temporary field and set it to resampler
    resampleFixed->SetTransform(negativeTrsf);
    resampleFixed->SetInput(m_BackwardImage);

    resampleFixed->SetSize(m_BackwardImage->GetLargestPossibleRegion().GetSize());
    resampleFixed->SetOutputOrigin(m_BackwardImage->GetOrigin());
    resampleFixed->SetOutputSpacing(m_BackwardImage->GetSpacing());
    resampleFixed->SetOutputDirection(m_BackwardImage->GetDirection());
    resampleFixed->SetDefaultPixelValue(0);
    resampleFixed->SetScaleIntensitiesWithJacobian(true);

    resampleFixed->Update();

    refImage = resampleFixed->GetOutput();
    refImage->DisconnectPipeline();
}

template <typename TInputImage>
void
DistortionCorrectionRegistrationMethod<TInputImage>
::InitializeBlocksOnImage(InputImageType *image)
{
    // Block creation should be moved here
    typedef anima::BlockMatchingInitializer<typename InputImageType::PixelType,
            TInputImage::ImageDimension> InitializerType;
    typename InitializerType::Pointer initPtr = InitializerType::New();

    initPtr->AddReferenceImage(image);
    initPtr->SetNumberOfThreads(this->GetNumberOfThreads());

    initPtr->SetPercentageKept(m_BlockPercentageKept);
    initPtr->SetBlockSize(m_BlockSize);
    initPtr->SetBlockSpacing(m_BlockSpacing);
    initPtr->SetScalarVarianceThreshold(m_BlockScalarVarianceThreshold);

    initPtr->SetRequestedRegion(image->GetLargestPossibleRegion());

    initPtr->SetComputeOuterDam(m_UseTransformationDam);
    initPtr->SetDamDistance(m_DamDistance);

    m_BlockRegions = initPtr->GetOutput();
    m_DamIndexes = initPtr->GetDamIndexes();

    m_Agregator->SetInputRegions(image, m_BlockRegions);

    std::cout << "Generated " << m_BlockRegions.size() << " blocks to process" << std::endl;

    m_BlockWeights.resize(m_BlockRegions.size());
    m_BaseTransformsPointers.resize(m_BlockRegions.size());
}

template <typename TInputImage>
void
DistortionCorrectionRegistrationMethod<TInputImage>
::PerformOneIteration(InputImageType *refImage, InputImageType *movingImage, SVFTransformType *addOn)
{
    this->InitializeBlocksOnImage(refImage);

    std::fill(m_BlockWeights.begin(),m_BlockWeights.end(),1.);
    // Do a threaded block match computation
    itk::MultiThreader::Pointer threader = itk::MultiThreader::New();

    ThreadedMatchData data;
    data.BlockMatch = this;
    data.fixedImage = refImage;
    data.movingImage = movingImage;

    itk::TimeProbe tmpTimeBackward;
    tmpTimeBackward.Start();

    unsigned int numThreads = std::min(this->GetNumberOfThreads(),(unsigned int)m_BlockWeights.size());

    threader->SetNumberOfThreads(numThreads);
    threader->SetSingleMethod(this->ThreadedMatching, &data);
    threader->SingleMethodExecute();

    tmpTimeBackward.Stop();
    std::cout << "Backward matching performed in " << tmpTimeBackward.GetTotal() << std::endl;

    // Do the agregation of solutions (positive then negative)
    m_Agregator->SetInputWeights(this->m_BlockWeights);
    m_Agregator->SetInputTransforms(this->m_BaseTransformsPointers);

    if (m_AgregatorType == Baloo)
        ((BalooSVFTransformAgregator<TInputImage::ImageDimension> *)m_Agregator)->SetDamIndexes(m_DamIndexes);
    else
        ((DenseSVFTransformAgregator<TInputImage::ImageDimension> *)m_Agregator)->SetDamIndexes(m_DamIndexes);

    typedef typename SVFTransformType::VectorFieldType VectorFieldType;
    SVFTransformType * tmpTrsf = (dynamic_cast <SVFTransformType *> (m_Agregator->GetOutput()));
    typename VectorFieldType::Pointer negativeSVF = const_cast <VectorFieldType *> (tmpTrsf->GetParametersAsVectorField());
    negativeSVF->DisconnectPipeline();

    this->InitializeBlocksOnImage(movingImage);

    std::fill(m_BlockWeights.begin(),m_BlockWeights.end(),1.);
    // Do a threaded block match computation
    threader = itk::MultiThreader::New();

    data.BlockMatch = this;
    data.fixedImage = movingImage;
    data.movingImage = refImage;

    itk::TimeProbe tmpTimeForward;
    tmpTimeForward.Start();

    threader->SetNumberOfThreads(numThreads);
    threader->SetSingleMethod(this->ThreadedMatching, &data);
    threader->SingleMethodExecute();

    tmpTimeForward.Stop();
    std::cout << "Forward matching performed in " << tmpTimeForward.GetTotal() << std::endl;

    // Do the agregation of solutions (positive then negative)
    m_Agregator->SetInputWeights(this->m_BlockWeights);
    m_Agregator->SetInputTransforms(this->m_BaseTransformsPointers);

    if (m_AgregatorType == Baloo)
        ((BalooSVFTransformAgregator<TInputImage::ImageDimension> *)m_Agregator)->SetDamIndexes(m_DamIndexes);
    else
        ((DenseSVFTransformAgregator<TInputImage::ImageDimension> *)m_Agregator)->SetDamIndexes(m_DamIndexes);

    tmpTrsf = (dynamic_cast <SVFTransformType *> (m_Agregator->GetOutput()));
    typename VectorFieldType::Pointer positiveSVF = const_cast <VectorFieldType *> (tmpTrsf->GetParametersAsVectorField());
    positiveSVF->DisconnectPipeline();

    typedef itk::MultiplyImageFilter <VectorFieldType,itk::Image <float,TInputImage::ImageDimension>,VectorFieldType> MultiplyFilterType;
    typedef itk::SubtractImageFilter <VectorFieldType,VectorFieldType,VectorFieldType> SubtractFilterType;

    typename SubtractFilterType::Pointer subFilter = SubtractFilterType::New();
    subFilter->SetInput1(positiveSVF);
    subFilter->SetInput2(negativeSVF);
    subFilter->SetNumberOfThreads(this->GetNumberOfThreads());
    subFilter->InPlaceOn();

    subFilter->Update();

    typename MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
    multiplyFilter->SetInput(subFilter->GetOutput());
    multiplyFilter->SetConstant(0.25);
    multiplyFilter->SetNumberOfThreads(this->GetNumberOfThreads());
    multiplyFilter->InPlaceOn();

    multiplyFilter->Update();

    positiveSVF = multiplyFilter->GetOutput();
    positiveSVF->DisconnectPipeline();

    addOn->SetParametersAsVectorField(positiveSVF);
}

/**
     * Threaded block match
     */
template <typename TInputImage>
ITK_THREAD_RETURN_TYPE
DistortionCorrectionRegistrationMethod<TInputImage>
::ThreadedMatching(void *arg)
{
    itk::MultiThreader::ThreadInfoStruct *threadArgs = (itk::MultiThreader::ThreadInfoStruct *)arg;
    // Set the registration procedure
    typename Self::ThreadedMatchData* data = (typename Self::ThreadedMatchData*)threadArgs->UserData;

    data->BlockMatch->BlockMatch(threadArgs->ThreadID, threadArgs->NumberOfThreads, data->fixedImage, data->movingImage);

    return NULL;
}

template <typename TInputImage>
void
DistortionCorrectionRegistrationMethod<TInputImage>
::GlobalParametersSetup()
{
    typedef anima::DirectionScaleSkewTransform <typename AgregatorType::ScalarType> BaseTransformType;
    typename BaseTransformType::Pointer tr;
    switch (m_TransformationKind)
    {
        case DirectionScaleSkew:
            tr = BaseTransformType::New();
            break;

        case DirectionScale:
            tr = anima::DirectionScaleTransform <typename AgregatorType::ScalarType>::New();
            break;

        case Direction:
        default:
            tr = anima::DirectionTransform <typename AgregatorType::ScalarType>::New();
            break;
    }

    tr->SetIdentity();
    m_InitialTransformParameters = tr->GetParameters();

    if (m_MetricKind == MeanSquares)
        m_MaximiseMetric = false;
    else
        m_MaximiseMetric = true;
}

template <typename TInputImage>
void
DistortionCorrectionRegistrationMethod<TInputImage>
::setUpParameters(InterpolatorPointer& interpolator,
                  MetricPointer&       metric,
                  OptimizerPointer&    optimizer,
                  unsigned startBlock,
                  unsigned endBlock)

{
    if (InputImageType::ImageDimension != 3)
        itkExceptionMacro("Only 3D transforms handled right now.");

    typedef anima::DirectionScaleSkewTransform <typename AgregatorType::ScalarType> BaseTransformType;
    typename BaseTransformType::HomogeneousMatrixType geometry;

    geometry.SetIdentity();
    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j < 3;++j)
            geometry(i,j) = m_BackwardImage->GetDirection()(i,j) * m_BackwardImage->GetSpacing()[j];

    for (unsigned i = startBlock; i < endBlock; i++)
    {
        typename BaseTransformType::Pointer tmpTr;

        switch (m_TransformationKind)
        {
            case DirectionScaleSkew:
                tmpTr = BaseTransformType::New();
                break;

            case DirectionScale:
                tmpTr = anima::DirectionScaleTransform <typename AgregatorType::ScalarType>::New();
                break;

            case Direction:
            default:
                tmpTr = anima::DirectionTransform <typename AgregatorType::ScalarType>::New();
                break;
        }

        tmpTr->SetIdentity();
        for (unsigned int j = 0;j < 3;++j)
            geometry(j,InputImageType::ImageDimension) = this->m_Agregator->GetInputOrigin(i)[j];

        tmpTr->SetUniqueDirection(m_TransformDirection);
        tmpTr->SetGeometry(geometry);

        m_BaseTransformsPointers[i] = tmpTr;
    }


    switch (m_MetricKind)
    {
        case MeanSquares:
        {
            typedef anima::FastMeanSquaresImageToImageMetric <InputImageType,InputImageType> LocalMetricType;

            typename LocalMetricType::Pointer tmpMetric = LocalMetricType::New();
            tmpMetric->SetScaleIntensities(true);
            metric = tmpMetric;
            break;
        }

        default:
        {
            typedef anima::FastCorrelationImageToImageMetric <InputImageType,InputImageType> LocalMetricType;
            
            typename LocalMetricType::Pointer tmpMetric = LocalMetricType::New();
            tmpMetric->SetScaleIntensities(true);
            tmpMetric->SetSquaredCorrelation(m_MetricKind == SquaredCorrelation);
            metric = tmpMetric;
            break;
        }
    }

    typedef anima::FasterLinearInterpolateImageFunction<InputImageType,double> LocalInterpolatorType;
    interpolator = LocalInterpolatorType::New();
    metric->SetInterpolator(interpolator);

    typedef anima::BobyqaOptimizer LocalOptimizerType;
    optimizer = LocalOptimizerType::New();
    LocalOptimizerType * tmpOpt = (LocalOptimizerType *)optimizer.GetPointer();
    tmpOpt->SetRhoBegin(m_SearchRadius);
    tmpOpt->SetRhoEnd(m_FinalRadius);

    tmpOpt->SetNumberSamplingPoints(m_BaseTransformsPointers[startBlock]->GetNumberOfParameters() + 2);
    // Maximum sampling of cost function ((Par+1)(Par+2))/2
    tmpOpt->SetMaximumIteration(m_OptimizerMaximumIterations);

    LocalOptimizerType::ScalesType tmpScales(m_BaseTransformsPointers[startBlock]->GetNumberOfParameters());
    LocalOptimizerType::ScalesType lowerBounds(m_BaseTransformsPointers[startBlock]->GetNumberOfParameters());
    LocalOptimizerType::ScalesType upperBounds(m_BaseTransformsPointers[startBlock]->GetNumberOfParameters());

    // Scale factor to ensure that max translations and skew can be reached
    // Based on the fact that non diagonal terms log is a = x * log(y) / (exp(y) - 1)
    // where y is the diagonal scaling factor, x the desired term
    double scaleFactor = 1.0;
    if ((m_TransformationKind != Direction)&&(m_ScaleMax > 0))
        scaleFactor = - m_ScaleMax / (std::exp(- m_ScaleMax) - 1.0);
    
    switch (m_TransformationKind)
    {
        case DirectionScaleSkew:
        {
            tmpScales[0] = m_SearchRadius / m_SearchScaleRadius;
            tmpScales[1] = m_SearchRadius / m_SearchSkewRadius;
            tmpScales[2] = m_SearchRadius / m_SearchSkewRadius;
            tmpScales[3] = 1.0;

            lowerBounds[0] = - m_ScaleMax;
            upperBounds[0] = m_ScaleMax;
            lowerBounds[1] = - m_SkewMax * scaleFactor;
            upperBounds[1] = m_SkewMax * scaleFactor;
            lowerBounds[2] = - m_SkewMax * scaleFactor;
            upperBounds[2] = m_SkewMax * scaleFactor;
            lowerBounds[3] = - m_TranslateMax * scaleFactor;
            upperBounds[3] = m_TranslateMax * scaleFactor;

            break;
        }

        case DirectionScale:
        {
            tmpScales[0] = m_SearchRadius / m_SearchScaleRadius;
            tmpScales[1] = 1.0;

            lowerBounds[0] = - m_ScaleMax;
            upperBounds[0] = m_ScaleMax;
            lowerBounds[1] = - m_TranslateMax * scaleFactor;
            upperBounds[1] = m_TranslateMax * scaleFactor;

            break;
        }


        case Direction:
        default:
        {
            tmpScales[0] = 1.0;
            lowerBounds[0] = - m_TranslateMax;
            upperBounds[0] = m_TranslateMax;

            break;
        }
    }

    tmpOpt->SetScales(tmpScales);
    tmpOpt->SetLowerBounds(lowerBounds);
    tmpOpt->SetUpperBounds(upperBounds);
    tmpOpt->SetMaximize(m_MaximiseMetric);
    metric->ComputeGradientOff();
}

template <typename TInputImage>
void
DistortionCorrectionRegistrationMethod<TInputImage>
::BlockMatch(unsigned threadId, unsigned NbThreads, InputImageType *fixedImage,
             InputImageType *movingImage)
{
    InterpolatorPointer interpolator;
    MetricPointer       metric;
    OptimizerPointer    optimizer;

    unsigned totalRegionSize = m_BlockRegions.size();
    unsigned step = totalRegionSize / NbThreads;

    unsigned startBlock = threadId*step;
    unsigned endBlock = (1+threadId)*step;

    if (threadId+1==NbThreads) // Ensure we got up to the last block
        endBlock = totalRegionSize;

    setUpParameters(interpolator, metric, optimizer, startBlock, endBlock);
    metric->SetFixedImage(fixedImage);
    metric->SetMovingImage(movingImage);

    if (GetDebug())
    {
        CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
        optimizer->AddObserver( itk::IterationEvent(), observer );
    }

    itkDebugMacro(<<"Threaded Block Match " << threadId << " / " << NbThreads << " starting @ " << startBlock << " up to " << endBlock );

    // Loop over the desired blocks
    typename BaseInputTransformType::Pointer realTransform;
    for (unsigned int block = startBlock; block < endBlock; ++block)
    {
        metric->SetFixedImageRegion(m_BlockRegions[block]);
        metric->SetTransform(m_BaseTransformsPointers[block]);
        metric->Initialize();

        optimizer->SetCostFunction(metric);
        optimizer->SetInitialPosition(m_InitialTransformParameters);

        if (m_MetricKind != MeanSquares)
            ((anima::FastCorrelationImageToImageMetric<TInputImage, TInputImage> *)metric.GetPointer())->PreComputeFixedValues();
        else
            ((anima::FastMeanSquaresImageToImageMetric<TInputImage, TInputImage> *)metric.GetPointer())->PreComputeFixedValues();

        try
        {
            optimizer->StartOptimization();
        }
        catch( itk::ExceptionObject & err )
        {
            std::cerr << "ExceptionObject caught Image Block Match! " << err << std::endl;
            m_BlockWeights[block] = 0;
            continue;
        }

        m_BaseTransformsPointers[block]->SetParameters( optimizer->GetCurrentPosition() );

        if (m_WeightedAgregation && (m_MetricKind != MeanSquares))
        {
            double val = ((anima::BobyqaOptimizer *)optimizer.GetPointer())->GetValue();
            m_BlockWeights[block] = (val >= 0) ? val : -val;
        }
        else
            m_BlockWeights[block] = 1;

        m_BlockWeights[block] *= this->ComputeDirectionFilterBlockWeightUpdate(fixedImage,m_BlockRegions[block]);
        m_BlockWeights[block] = std::pow(m_BlockWeights[block],1/2.0);

        itkDebugMacro(<<  m_BaseTransformsPointers[block]->GetParameters() << "  " << m_BlockWeights[block] << std::endl);
    }
}

template <typename TInputImage>
double
DistortionCorrectionRegistrationMethod<TInputImage>
::ComputeDirectionFilterBlockWeightUpdate(InputImageType *fixedImage, const ImageRegionType &region)
{
    std::vector <double> localGradient(TInputImage::ImageDimension,0);
    itk::ImageRegionConstIterator <InputImageType> fixedItr(fixedImage,region);
    typename ImageRegionType::IndexType currentIndex, modifiedIndex;

    typename InputImageType::DirectionType orientationMatrix = fixedImage->GetDirection();
    typename InputImageType::SpacingType imageSpacing = fixedImage->GetSpacing();
    typename InputImageType::SizeType imageSize = fixedImage->GetLargestPossibleRegion().GetSize();

    std::vector <double> correctionDirection(TInputImage::ImageDimension);
    for (unsigned int i = 0;i < TInputImage::ImageDimension;++i)
        correctionDirection[i] = fixedImage->GetDirection()(i,m_TransformDirection);

    vnl_matrix <double> meanStructureTensor(TInputImage::ImageDimension,TInputImage::ImageDimension);
    meanStructureTensor.fill(0);

    while (!fixedItr.IsAtEnd())
    {
        currentIndex = fixedItr.GetIndex();
        std::fill(localGradient.begin(),localGradient.end(),0);

        for (unsigned int i = 0;i < TInputImage::ImageDimension;++i)
        {
            modifiedIndex = currentIndex;
            modifiedIndex[i] = std::max(0,(int)(currentIndex[i] - 1));
            double previousValue = fixedImage->GetPixel(modifiedIndex);
            modifiedIndex[i] = std::min((int)(imageSize[i] - 1),(int)(currentIndex[i] + 1));
            double postValue = fixedImage->GetPixel(modifiedIndex);

            for (unsigned int j = 0;j < TInputImage::ImageDimension;++j)
                localGradient[j] += (postValue - previousValue) * orientationMatrix(j,i) / (2.0 * imageSpacing[i]);
        }

        for (unsigned int i = 0;i < TInputImage::ImageDimension;++i)
            for (unsigned int j = i;j < TInputImage::ImageDimension;++j)
            {
                meanStructureTensor(i,j) += localGradient[i] * localGradient[j];
                if (j != i)
                    meanStructureTensor(j,i) = meanStructureTensor(i,j);
            }

        ++fixedItr;
    }

    meanStructureTensor /= region.GetNumberOfPixels();

    itk::SymmetricEigenAnalysis < vnl_matrix <double>, vnl_diag_matrix<double>, vnl_matrix <double> > eigenComputer(TInputImage::ImageDimension);
    vnl_matrix <double> eVec(TInputImage::ImageDimension,TInputImage::ImageDimension);
    vnl_diag_matrix <double> eVals(TInputImage::ImageDimension);

    eigenComputer.ComputeEigenValuesAndVectors(meanStructureTensor, eVals, eVec);
    double linearCoef = (eVals[TInputImage::ImageDimension - 1] - eVals[TInputImage::ImageDimension - 2]) / eVals[TInputImage::ImageDimension - 1];

    double scalarProduct = 0;
    for (unsigned int i = 0;i < TInputImage::ImageDimension;++i)
        scalarProduct += eVec[TInputImage::ImageDimension - 1][i] * correctionDirection[i];

    return linearCoef * std::abs(scalarProduct);
    //return std::abs(scalarProduct);
}

/**
     * PrintSelf
     */
template <typename TInputImage>
void
DistortionCorrectionRegistrationMethod<TInputImage>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
    Superclass::PrintSelf( os, indent );
    os << indent << "Backward Image: " << m_BackwardImage.GetPointer() << std::endl;
    os << indent << "Forward Image: " << m_ForwardImage.GetPointer() << std::endl;

    os << indent << "Search Radius: " << m_SearchRadius << std::endl;
    os << indent << "Final Radius: " << m_FinalRadius << std::endl;
    os << indent << "Step Size : " << m_StepSize << std::endl;

    os << indent << "Weighted Agregation : " << m_WeightedAgregation << std::endl;

    os << indent << "Maximum Iterations: " << m_MaximumIterations << std::endl;

    os << indent << "Optimizer Maximum Iterations: " << m_OptimizerMaximumIterations << std::endl;

    os << indent << "Initial Transform Parameters: " << m_InitialTransformParameters << std::endl;
    os << indent << "Last    Transform Parameters: " << m_LastTransformParameters << std::endl;
}

/*
     * Generate Data
     */
template <typename TInputImage>
void
DistortionCorrectionRegistrationMethod<TInputImage>
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
    catch( itk::ExceptionObject& err )
    {
        m_LastTransformParameters = empty;

        // pass exception to caller
        throw err;
    }

    this->StartOptimization();
}

/**
     *  Get Output
     */
template <typename TInputImage>
typename DistortionCorrectionRegistrationMethod<TInputImage>::TransformOutputType *
DistortionCorrectionRegistrationMethod<TInputImage>
::GetOutput()
{
    return static_cast< TransformOutputType * >( this->ProcessObject::GetOutput(0) );
}

template <typename TInputImage>
itk::DataObject::Pointer
DistortionCorrectionRegistrationMethod<TInputImage>
::MakeOutput(DataObjectPointerArraySizeType output)
{
    switch (output)
    {
        case 0:
            return static_cast<itk::DataObject*>(TransformOutputType::New().GetPointer());
            break;
        default:
            itkExceptionMacro("MakeOutput request for an output number larger than the expected number of outputs");
            return 0;
    }
}

template <typename TInputImage>
void
DistortionCorrectionRegistrationMethod<TInputImage>
::SetBackwardImage( InputImageType * backwardImage )
{
    itkDebugMacro("setting backward Image to " << backwardImage );

    if (this->m_BackwardImage.GetPointer() != backwardImage )
    {
        this->m_BackwardImage = backwardImage;
        this->m_BackwardImage->DisconnectPipeline();

        this->Modified();
    }
}

template <typename TInputImage>
void
DistortionCorrectionRegistrationMethod<TInputImage>
::SetForwardImage( InputImageType * forwardImage )
{
    itkDebugMacro("setting forward Image to " << forwardImage );

    if (this->m_ForwardImage.GetPointer() != forwardImage )
    {
        this->m_ForwardImage = forwardImage;
        this->m_ForwardImage->DisconnectPipeline();

        this->Modified();
    }
}

template <typename TInputImage>
const std::string
DistortionCorrectionRegistrationMethod<TInputImage>
::GetStopConditionDescription() const
{
    return m_StopConditionDescription.str();
}

} // end namespace anima
