#pragma once
#include <animaMCMBlockMatcher.h>

/* Similarity measures */
#include <animaMCMBasicMeanSquaresImageToImageMetric.h>
#include <animaMCMMeanSquaresImageToImageMetric.h>
#include <animaMCMCorrelationImageToImageMetric.h>
#include <animaBaseOrientedModelImageToImageMetric.h>

#include <animaMCMBlockMatchInitializer.h>

namespace anima
{

template <typename TInputImageType>
MCMBlockMatcher<TInputImageType>
::MCMBlockMatcher()
{
    m_SimilarityType = MCMMeanSquares;
}

template <typename TInputImageType>
bool
MCMBlockMatcher<TInputImageType>
::GetMaximizedMetric()
{
    return false;
}

template <typename TInputImageType>
typename MCMBlockMatcher<TInputImageType>::MCMInterpolatorType *
MCMBlockMatcher<TInputImageType>
::CreateInterpolator()
{
    typename MCMInterpolatorType::Pointer interpolator = MCMInterpolatorType::New();
    interpolator->SetReferenceOutputModel(this->GetReferenceImage()->GetDescriptionModel());

    interpolator->Register();
    return interpolator;
}

template <typename TInputImageType>
void
MCMBlockMatcher<TInputImageType>
::InitializeBlocks()
{
    typedef typename TInputImageType::IOPixelType InputPixelType;
    typedef typename anima::MCMBlockMatchingInitializer<InputPixelType,TInputImageType::ImageDimension> InitializerType;
    typedef typename InitializerType::Pointer InitializerPointer;

    InitializerPointer initPtr = InitializerType::New();
    initPtr->AddReferenceImage(this->GetReferenceImage());

    if (this->GetNumberOfThreads() != 0)
        initPtr->SetNumberOfThreads(this->GetNumberOfThreads());

    initPtr->SetPercentageKept(this->GetBlockPercentageKept());
    initPtr->SetBlockSize(this->GetBlockSize());
    initPtr->SetBlockSpacing(this->GetBlockSpacing());
    initPtr->SetOrientedModelVarianceThreshold(this->GetBlockVarianceThreshold());

    initPtr->SetRequestedRegion(this->GetReferenceImage()->GetLargestPossibleRegion());

    initPtr->SetComputeOuterDam(this->GetUseTransformationDam());
    initPtr->SetDamDistance(this->GetDamDistance());

    this->SetBlockRegions(initPtr->GetOutput());
    this->SetBlockPositions(initPtr->GetOutputPositions());
    this->SetBlockDamWeights(initPtr->GetBlockDamWeights());

    std::cout << "Generated " << this->GetBlockRegions().size() << " blocks..." << std::endl;

    this->GetBlockTransformPointers().resize(this->GetBlockRegions().size());
    std::vector <double> newBlockWeights(this->GetBlockRegions().size(),0);
    this->SetBlockWeights(newBlockWeights);
    for (unsigned int i = 0;i < this->GetBlockRegions().size();++i)
        this->GetBlockTransformPointer(i) = this->GetNewBlockTransform(this->GetBlockPositions()[i]);
}

template <typename TInputImageType>
typename MCMBlockMatcher<TInputImageType>::MetricPointer
MCMBlockMatcher<TInputImageType>
::SetupMetric()
{
    MetricPointer metric;

    switch(m_SimilarityType)
    {
        case MCMBasicMeanSquares:
        case MCMOneToOneBasicMeanSquares:
        {
            typedef anima::MCMBasicMeanSquaresImageToImageMetric<typename InputImageType::IOPixelType, typename InputImageType::IOPixelType,
                                                                 InputImageType::ImageDimension > MetricType;

            typename MetricType::Pointer tmpMetric = MetricType::New();
            tmpMetric->SetOneToOneMapping(m_SimilarityType == MCMOneToOneBasicMeanSquares);
            metric = tmpMetric;
            break;
        }

        case MCMCorrelation:
        {
            typedef anima::MCMCorrelationImageToImageMetric<typename InputImageType::IOPixelType, typename InputImageType::IOPixelType,
                                                            InputImageType::ImageDimension > MetricType;

            typename MetricType::Pointer tmpMetric = MetricType::New();
            tmpMetric->SetBValues(m_BValues);
            tmpMetric->SetGradientDirections(m_GradientDirections);
            if ((m_BValues.size() > 0)&&(m_GradientDirections.size() > 0))
                tmpMetric->SetForceApproximation(true);

            metric = tmpMetric;
            break;
        }

        case MCMMeanSquares:
        default:
        {
            typedef anima::MCMMeanSquaresImageToImageMetric<typename InputImageType::IOPixelType, typename InputImageType::IOPixelType,
                                                            InputImageType::ImageDimension > MetricType;

            typename MetricType::Pointer tmpMetric = MetricType::New();
            tmpMetric->SetBValues(m_BValues);
            tmpMetric->SetGradientDirections(m_GradientDirections);
            if ((m_BValues.size() > 0)&&(m_GradientDirections.size() > 0))
                tmpMetric->SetForceApproximation(true);

            metric = tmpMetric;
            break;
        }
    }

    typedef anima::BaseOrientedModelImageToImageMetric <InputImageType,InputImageType> BaseMetricType;
    BaseMetricType *baseMetric = dynamic_cast <BaseMetricType *> (metric.GetPointer());

    if (this->GetBlockTransformType() == Superclass::Translation)
        baseMetric->SetModelRotation(BaseMetricType::NONE);
    else
        baseMetric->SetModelRotation((typename BaseMetricType::ModelReorientationType)m_ModelRotationType);

    typename MCMInterpolatorType::Pointer interpolator = this->CreateInterpolator();

    baseMetric->SetInterpolator(interpolator);

    baseMetric->SetFixedImage(this->GetReferenceImage());
    baseMetric->SetMovingImage(this->GetMovingImage());

    return metric;
}

template <typename TInputImageType>
double
MCMBlockMatcher<TInputImageType>
::ComputeBlockWeight(double val, unsigned int block)
{    
    return 1.0;
}

template <typename TInputImageType>
void
MCMBlockMatcher<TInputImageType>
::BlockMatchingSetup(MetricPointer &metric, unsigned int block)
{
    // For transform related init
    this->Superclass::BlockMatchingSetup(metric,block);

    // Metric specific init
    typedef anima::BaseOrientedModelImageToImageMetric <InputImageType,InputImageType> InternalMetricType;
    InternalMetricType *tmpMetric = dynamic_cast <InternalMetricType *> (metric.GetPointer());
    tmpMetric->SetFixedImageRegion(this->GetBlockRegion(block));
    tmpMetric->SetTransform(this->GetBlockTransformPointer(block));
    tmpMetric->Initialize();

    typedef anima::MCMBasicMeanSquaresImageToImageMetric <typename InputImageType::IOPixelType,
                                                          typename InputImageType::IOPixelType,
                                                          InputImageType::ImageDimension > MCMBasicMeanSquaresMetricType;

    typedef anima::MCMMeanSquaresImageToImageMetric <typename InputImageType::IOPixelType,
                                                     typename InputImageType::IOPixelType,
                                                     InputImageType::ImageDimension > MCMMeanSquaresMetricType;

    typedef anima::MCMCorrelationImageToImageMetric <typename InputImageType::IOPixelType,
                                                     typename InputImageType::IOPixelType,
                                                     InputImageType::ImageDimension > MCMCorrelationMetricType;

    switch(m_SimilarityType)
    {
        case MCMMeanSquares:
            ((MCMMeanSquaresMetricType *)metric.GetPointer())->PreComputeFixedValues();
            break;
        case MCMCorrelation:
            ((MCMCorrelationMetricType *)metric.GetPointer())->PreComputeFixedValues();
            break;
        case MCMBasicMeanSquares:
        case MCMOneToOneBasicMeanSquares:
        default:
            ((MCMBasicMeanSquaresMetricType *)metric.GetPointer())->PreComputeFixedValues();
            break;
    }
}

} // end namespace anima
