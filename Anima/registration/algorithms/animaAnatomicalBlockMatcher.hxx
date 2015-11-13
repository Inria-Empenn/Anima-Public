#pragma once
#include <animaAnatomicalBlockMatcher.h>

/* Similarity measures */
#include <animaFastCorrelationImageToImageMetric.h>
#include <itkMeanSquaresImageToImageMetric.h>
#include <itkImageToImageMetric.h>

#include <animaFasterLinearInterpolateImageFunction.h>

namespace anima
{

template <typename TInputImageType>
AnatomicalBlockMatcher<TInputImageType>
::AnatomicalBlockMatcher()
{
    m_SimilarityType = SquaredCorrelation;
}

template <typename TInputImageType>
bool
AnatomicalBlockMatcher<TInputImageType>
::GetMaximizedMetric()
{
    if (m_SimilarityType == MeanSquares)
        return false;

    return true;
}

template <typename TInputImageType>
typename AnatomicalBlockMatcher<TInputImageType>::MetricPointer
AnatomicalBlockMatcher<TInputImageType>
::SetupMetric()
{
    MetricPointer metric;

    switch(m_SimilarityType)
    {
        case Correlation:
        case SquaredCorrelation:
        {
            typedef anima::FastCorrelationImageToImageMetric <InputImageType,InputImageType> LocalMetricType;

            typename LocalMetricType::Pointer tmpMetric = LocalMetricType::New();
            tmpMetric->SetSquaredCorrelation(m_SimilarityType == SquaredCorrelation);

            metric = tmpMetric;
            break;
        }

        case MeanSquares:
        default:
        {
            typedef itk::MeanSquaresImageToImageMetric <InputImageType,InputImageType> LocalMetricType;

            metric = LocalMetricType::New();
            break;
        }
    }

    typedef itk::ImageToImageMetric <InputImageType,InputImageType> BaseMetricType;
    BaseMetricType *baseMetric = dynamic_cast <BaseMetricType *> (metric.GetPointer());

    typedef anima::FasterLinearInterpolateImageFunction<InputImageType,double> LocalInterpolatorType;
    typename LocalInterpolatorType::Pointer interpolator = LocalInterpolatorType::New();

    baseMetric->SetInterpolator(interpolator);
    baseMetric->ComputeGradientOff();

    baseMetric->SetFixedImage(this->GetReferenceImage());
    baseMetric->SetMovingImage(this->GetMovingImage());
    interpolator->SetInputImage(this->GetMovingImage());

    return metric;
}

template <typename TInputImageType>
double
AnatomicalBlockMatcher<TInputImageType>
::ComputeBlockWeight(double val, unsigned int block)
{
    switch (m_SimilarityType)
    {
        case MeanSquares:
            return 1;

        case Correlation:
            return (val + 1) / 2.0;

        case SquaredCorrelation:
        default:
            return val;
    }
}

template <typename TInputImageType>
void
AnatomicalBlockMatcher<TInputImageType>
::BlockMatchingSetup(MetricPointer &metric, unsigned int block)
{
    // For transform related init
    this->Superclass::BlockMatchingSetup(metric,block);

    // Metric specific init
    typedef itk::ImageToImageMetric <InputImageType, InputImageType> InternalMetricType;
    InternalMetricType *tmpMetric = dynamic_cast <InternalMetricType *> (metric.GetPointer());
    tmpMetric->SetFixedImageRegion(this->GetBlockRegion(block));
    tmpMetric->SetTransform(this->GetBlockTransformPointer(block));
    tmpMetric->Initialize();
    if (m_SimilarityType != MeanSquares)
        ((anima::FastCorrelationImageToImageMetric<InputImageType, InputImageType> *)metric.GetPointer())->PreComputeFixedValues();
}

} // end namespace anima
