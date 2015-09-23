#pragma once
#include <animaTensorBlockMatcher.h>

/* Similarity measures */
#include <animaTensorCorrelationImageToImageMetric.h>
#include <animaTensorGeneralizedCorrelationImageToImageMetric.h>
#include <animaTensorMeanSquaresImageToImageMetric.h>
#include <animaBaseTensorImageToImageMetric.h>

#include <animaVectorModelLinearInterpolateImageFunction.h>

namespace anima
{

template <typename TInputImageType>
TensorBlockMatcher<TInputImageType>
::TensorBlockMatcher()
{
    m_SimilarityType = TensorOrientedGeneralizedCorrelation;
}

template <typename TInputImageType>
bool
TensorBlockMatcher<TInputImageType>
::GetMaximizedMetric()
{
    if (m_SimilarityType == TensorMeanSquares)
        return false;

    return true;
}

template <typename TInputImageType>
typename TensorBlockMatcher<TInputImageType>::MetricPointer
TensorBlockMatcher<TInputImageType>
::SetupMetric()
{
    MetricPointer metric;

    switch(m_SimilarityType)
    {
        case TensorMeanSquares:
        {
            typedef anima::TensorMeanSquaresImageToImageMetric<typename InputImageType::IOPixelType,
                                                             typename InputImageType::IOPixelType, InputImageType::ImageDimension > MetricType;

            metric = MetricType::New();
            break;
        }

        case TensorGeneralizedCorrelation:
        case TensorOrientedGeneralizedCorrelation:
        {
            typedef anima::TensorGeneralizedCorrelationImageToImageMetric<typename InputImageType::IOPixelType,
                                                                        typename InputImageType::IOPixelType, InputImageType::ImageDimension > MetricType;

            typename MetricType::Pointer tmpMetric = MetricType::New();

            tmpMetric->SetOrientationPenalty(m_SimilarityType == TensorOrientedGeneralizedCorrelation);
            tmpMetric->SetVarianceThreshold(this->GetBlockVarianceThreshold());

            metric = tmpMetric;
            break;
        }

        case TensorCorrelation:
        default:
        {
            typedef anima::TensorCorrelationImageToImageMetric<typename InputImageType::IOPixelType,
                                                             typename InputImageType::IOPixelType, InputImageType::ImageDimension > MetricType;

            metric = MetricType::New();
            break;
        }
    }

    typedef anima::BaseTensorImageToImageMetric <InputImageType,InputImageType> BaseMetricType;
    BaseMetricType *baseMetric = dynamic_cast <BaseMetricType *> (metric.GetPointer());
    baseMetric->SetRotateTensors(this->GetBlockTransformType() != Superclass::Translation);

    typedef anima::VectorModelLinearInterpolateImageFunction<InputImageType,double> LocalInterpolatorType;
    typename LocalInterpolatorType::Pointer interpolator = LocalInterpolatorType::New();

    baseMetric->SetInterpolator(interpolator);

    baseMetric->SetFixedImage(this->GetReferenceImage());
    baseMetric->SetMovingImage(this->GetMovingImage());
    interpolator->SetInputImage(this->GetMovingImage());

    return metric;
}

template <typename TInputImageType>
double
TensorBlockMatcher<TInputImageType>
::ComputeBlockWeight(double val, unsigned int block)
{
    switch (m_SimilarityType)
    {
        case TensorMeanSquares:
            return 1;

        default:
            return val;
    }
}

template <typename TInputImageType>
void
TensorBlockMatcher<TInputImageType>
::BlockMatchingSetup(MetricPointer &metric, unsigned int block)
{
    // For transform related init
    this->Superclass::BlockMatchingSetup(metric,block);

    // Metric specific init
    typedef anima::BaseTensorImageToImageMetric <InputImageType,InputImageType> InternalMetricType;
    InternalMetricType *tmpMetric = dynamic_cast <InternalMetricType *> (metric.GetPointer());
    tmpMetric->SetFixedImageRegion(this->GetBlockRegion(block));
    tmpMetric->SetTransform(this->GetBlockTransformPointer(block));
    tmpMetric->Initialize();

    typedef anima::TensorCorrelationImageToImageMetric <typename InputImageType::IOPixelType,
                                                      typename InputImageType::IOPixelType,
                                                      InputImageType::ImageDimension > TensorCorrelationMetricType;


    typedef anima::TensorGeneralizedCorrelationImageToImageMetric <typename InputImageType::IOPixelType,
                                                                 typename InputImageType::IOPixelType,
                                                                 InputImageType::ImageDimension > TensorGeneralizedMetricType;


    if ((m_SimilarityType == TensorGeneralizedCorrelation)||(m_SimilarityType == TensorOrientedGeneralizedCorrelation))
        ((TensorGeneralizedMetricType *)metric.GetPointer())->PreComputeFixedValues();

    if (m_SimilarityType == TensorCorrelation)
        ((TensorCorrelationMetricType *)metric.GetPointer())->PreComputeFixedValues();
}

} // end namespace anima
