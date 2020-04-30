#pragma once
#include "animaFastMeanSquaresImageToImageMetric.h"

#include <itkImageRegionConstIteratorWithIndex.h>

namespace anima
{

template <class TFixedImage, class TMovingImage>
FastMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::FastMeanSquaresImageToImageMetric()
{
    m_ScaleIntensities = false;
    m_DefaultBackgroundValue = 0.0;
}

template <class TFixedImage, class TMovingImage>
typename FastMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>::MeasureType
FastMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::GetValue( const TransformParametersType & parameters ) const
{
    FixedImageConstPointer fixedImage = this->m_FixedImage;

    if( !fixedImage )
    {
        itkExceptionMacro( << "Fixed image has not been assigned" );
    }

    if (this->m_NumberOfPixelsCounted == 0)
        return 0;

    MeasureType measure = 0;
    this->SetTransformParameters( parameters );

    OutputPointType transformedPoint;
    ContinuousIndexType transformedIndex;
    RealType movingValue;

    for (unsigned int i = 0;i < this->m_NumberOfPixelsCounted;++i)
    {
        transformedPoint = this->m_Transform->TransformPoint( m_FixedImagePoints[i] );
        this->m_Interpolator->GetInputImage()->TransformPhysicalPointToContinuousIndex(transformedPoint,transformedIndex);

        movingValue = m_DefaultBackgroundValue;

        if( this->m_Interpolator->IsInsideBuffer( transformedIndex ) )
        {
            movingValue = this->m_Interpolator->EvaluateAtContinuousIndex( transformedIndex );

            if (m_ScaleIntensities)
            {
                typedef itk::MatrixOffsetTransformBase <typename TransformType::ScalarType,
                        TFixedImage::ImageDimension, TFixedImage::ImageDimension> BaseTransformType;
                BaseTransformType *currentTrsf = dynamic_cast<BaseTransformType *> (this->m_Transform.GetPointer());

                double factor = vnl_determinant(currentTrsf->GetMatrix().GetVnlMatrix());
                movingValue *= factor;
            }
        }

        measure += (movingValue - m_FixedImageValues[i]) * (movingValue - m_FixedImageValues[i]);
    }

    measure /= this->m_NumberOfPixelsCounted;

    return measure;
}

template < class TFixedImage, class TMovingImage>
void
FastMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::GetDerivative( const TransformParametersType & parameters,
                 DerivativeType & derivative ) const
{
    itkExceptionMacro("Derivative not implemented yet...");
}

template <class TFixedImage, class TMovingImage>
void
FastMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::GetValueAndDerivative(const TransformParametersType & parameters,
                        MeasureType & value, DerivativeType  & derivative) const
{
    itkExceptionMacro("Derivative not implemented yet...");
}

template <class TFixedImage, class TMovingImage>
void
FastMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::PreComputeFixedValues()
{
    FixedImageConstPointer fixedImage = this->m_FixedImage;

    if( !fixedImage )
    {
        itkExceptionMacro( << "Fixed image has not been assigned" );
    }

    typedef itk::ImageRegionConstIteratorWithIndex<FixedImageType> FixedIteratorType;

    FixedIteratorType ti( fixedImage, this->GetFixedImageRegion() );
    typename FixedImageType::IndexType index;

    this->m_NumberOfPixelsCounted = this->GetFixedImageRegion().GetSize()[0];
    for (unsigned int i = 1;i < TFixedImage::GetImageDimension();++i)
        this->m_NumberOfPixelsCounted *= this->GetFixedImageRegion().GetSize()[i];

    m_FixedImagePoints.resize(this->m_NumberOfPixelsCounted);
    m_FixedImageValues.resize(this->m_NumberOfPixelsCounted);

    InputPointType inputPoint;

    unsigned int pos = 0;
    RealType fixedValue;

    while(!ti.IsAtEnd())
    {
        index = ti.GetIndex();
        fixedImage->TransformIndexToPhysicalPoint( index, inputPoint );

        m_FixedImagePoints[pos] = inputPoint;
        fixedValue = ti.Value();
        m_FixedImageValues[pos] = fixedValue;

        ++ti;
        ++pos;
    }
}

} // end namespace anima
