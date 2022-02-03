#pragma once
#include "animaFastCorrelationImageToImageMetric.h"

#include <itkImageRegionConstIteratorWithIndex.h>

namespace anima
{

template <class TFixedImage, class TMovingImage>
FastCorrelationImageToImageMetric<TFixedImage,TMovingImage>
::FastCorrelationImageToImageMetric()
{
    m_SumFixed = 0;
    m_VarFixed = 0;
    m_DefaultBackgroundValue = 0.0;
    m_SquaredCorrelation = true;
    m_ScaleIntensities = false;
    m_FixedImagePoints.clear();
    m_FixedImageValues.clear();
}

template <class TFixedImage, class TMovingImage>
typename FastCorrelationImageToImageMetric<TFixedImage,TMovingImage>::MeasureType
FastCorrelationImageToImageMetric<TFixedImage,TMovingImage>
::GetValue( const TransformParametersType & parameters ) const
{
    FixedImageConstPointer fixedImage = this->m_FixedImage;

    if (!fixedImage)
        itkExceptionMacro( << "Fixed image has not been assigned" );

    if ( this->m_NumberOfPixelsCounted == 0 )
        return 0;

    MeasureType measure;
    this->SetTransformParameters( parameters );

    typedef typename itk::NumericTraits< MeasureType >::AccumulateType AccumulateType;

    AccumulateType smm = itk::NumericTraits< AccumulateType >::Zero;
    AccumulateType sfm = itk::NumericTraits< AccumulateType >::Zero;
    AccumulateType sm  = itk::NumericTraits< AccumulateType >::Zero;

    OutputPointType transformedPoint;
    ContinuousIndexType transformedIndex;
    RealType movingValue;

    for (unsigned int i = 0;i < this->m_NumberOfPixelsCounted;++i)
    {
        transformedPoint = this->m_Transform->TransformPoint(m_FixedImagePoints[i]);
        this->m_Interpolator->GetInputImage()->TransformPhysicalPointToContinuousIndex(transformedPoint,transformedIndex);

        movingValue = m_DefaultBackgroundValue;
        if (this->m_Interpolator->IsInsideBuffer(transformedIndex))
            movingValue = this->m_Interpolator->EvaluateAtContinuousIndex(transformedIndex);

        if (movingValue != 0.0)
        {
            if (m_ScaleIntensities)
            {
                typedef itk::MatrixOffsetTransformBase <typename TransformType::ScalarType,
                        TFixedImage::ImageDimension, TFixedImage::ImageDimension> BaseTransformType;
                BaseTransformType *currentTrsf = dynamic_cast<BaseTransformType *> (this->m_Transform.GetPointer());

                double factor = vnl_determinant(currentTrsf->GetMatrix().GetVnlMatrix());
                movingValue *= factor;
            }

            smm += movingValue * movingValue;
            sfm += m_FixedImageValues[i] * movingValue;
            sm += movingValue;
        }
    }

    RealType movingVariance = smm - sm * sm / this->m_NumberOfPixelsCounted;
    RealType covData = sfm - m_SumFixed * sm / this->m_NumberOfPixelsCounted;
    RealType multVars = m_VarFixed * movingVariance;

    if (this->m_NumberOfPixelsCounted > 1 && multVars > 1.0e-16)
    {
        if (m_SquaredCorrelation)
            measure = covData * covData / multVars;
        else
            measure = std::max(0.0,covData / sqrt(multVars));
    }
    else
    {
        if (m_SquaredCorrelation)
            measure = itk::NumericTraits< MeasureType >::Zero;
        else
            measure = - 1.0;
    }

    return measure;
}

template < class TFixedImage, class TMovingImage>
void
FastCorrelationImageToImageMetric<TFixedImage,TMovingImage>
::PreComputeFixedValues()
{
    FixedImageConstPointer fixedImage = this->m_FixedImage;

    if( !fixedImage )
    {
        itkExceptionMacro( << "Fixed image has not been assigned" );
    }

    m_SumFixed = 0;
    m_VarFixed = 0;
    RealType sumSquared = 0;

    typedef itk::ImageRegionConstIteratorWithIndex<FixedImageType> FixedIteratorType;

    FixedIteratorType ti( fixedImage, this->GetFixedImageRegion() );
    typename FixedImageType::IndexType index;

    this->m_NumberOfPixelsCounted = this->GetFixedImageRegion().GetNumberOfPixels();

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

        sumSquared += fixedValue * fixedValue;
        m_SumFixed += fixedValue;

        ++ti;
        ++pos;
    }

    m_VarFixed = sumSquared - m_SumFixed * m_SumFixed / this->m_NumberOfPixelsCounted;
}

template < class TFixedImage, class TMovingImage>
void
FastCorrelationImageToImageMetric<TFixedImage,TMovingImage>
::GetDerivative( const TransformParametersType & parameters,
                DerivativeType & derivative ) const
{
    itkExceptionMacro("Derivative not implemented yet...");
}

template <class TFixedImage, class TMovingImage>
void
FastCorrelationImageToImageMetric<TFixedImage,TMovingImage>
::GetValueAndDerivative(const TransformParametersType & parameters,
                        MeasureType & value, DerivativeType  & derivative) const
{
    itkExceptionMacro("Derivative not implemented yet...");
}

template < class TFixedImage, class TMovingImage>
void
FastCorrelationImageToImageMetric<TFixedImage,TMovingImage>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
    os << indent << m_SumFixed << " " << m_VarFixed << std::endl;
}

} // end of namespace anima
