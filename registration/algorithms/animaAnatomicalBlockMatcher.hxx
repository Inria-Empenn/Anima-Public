#pragma once
#include <animaAnatomicalBlockMatcher.h>

/* Transforms */
#include <itkTranslationTransform.h>
#include <animaLogRigid3DTransform.h>
#include <animaSplitAffine3DTransform.h>

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
    m_BlockTransformType = Translation;
    m_SimilarityType = SquaredCorrelation;
}

template <typename TInputImageType>
typename AnatomicalBlockMatcher<TInputImageType>::AgregatorType::TRANSFORM_TYPE
AnatomicalBlockMatcher<TInputImageType>
::GetAgregatorInputTransformType()
{
    switch (m_BlockTransformType)
    {
        case Translation:
            return AgregatorType::TRANSLATION;

        case Rigid:
            return AgregatorType::RIGID;

        case Affine:
        default:
            return AgregatorType::AFFINE;
    }
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
AnatomicalBlockMatcher<TInputImageType>::MetricPointer
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
            tmpMetric->SetSquaredCorrelation(m_MetricKind == SquaredCorrelation);
            tmpMetric->SetVarianceThreshold(this->GetBlockVarianceThreshold());

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
    BaseMetricType *baseMetric = dynamic_cast <BaseMetricType *> (metric);

    typedef anima::FasterLinearInterpolateImageFunction<InputImageType,double> LocalInterpolatorType;
    interpolator = LocalInterpolatorType::New();

    baseMetric->SetInterpolator(interpolator);
    baseMetric->ComputeGradientOff();

    interpolator->SetInputImage(m_MovingImage);
}

template <typename TInputImageType>
double
AnatomicalBlockMatcher<TInputImageType>
::ComputeBlockWeight(double val)
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
typename AnatomicalBlockMatcher<TInputImageType>::BaseInputTransformPointer
AnatomicalBlockMatcher<TInputImageType>
::GetNewBlockTransform()
{
    BaseInputTransformPointer outputValue;

    switch(m_BlockTransformType)
    {
        case Translation:
        {
            typedef itk::TranslationTransform <double, InputImageType::ImageDimension> itkTransformType;
            typename itkTransformType::Pointer tr = itkTransformType::New();
            tr->SetIdentity();
            outputValue = tr;

            break;
        }

        case Rigid:
        {
            if (InputImageType::ImageDimension == 3)
            {
                typedef anima::LogRigid3DTransform<double> itkTransformType;
                typename itkTransformType::Pointer tr = itkTransformType::New();
                tr->SetIdentity();

                outputValue = tr;
            }
            else
                std::cerr << "Only Rigid 3D transforms handled right now." << std::endl;
            break;
        }

        case Affine:
        default:
        {
            typedef anima::SplitAffine3DTransform <double> itkTransformType;
            typename itkTransformType::Pointer tr = itkTransformType::New();
            tr->SetIdentity();

            outputValue = tr;

            break;
        }
    }

    return outputValue;
}

template <typename TInputImageType>
void
AnatomicalBlockMatcher<TInputImageType>
::BlockMatchingSetup(MetricPointer &metric, unsigned int block)
{

}

template <typename TInputImageType>
void
AnatomicalBlockMatcher<TInputImageType>
::TransformDependantOptimizerSetup(OptimizerPointer &optimizer)
{

}

} // end namespace anima
