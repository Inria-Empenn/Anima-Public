#pragma once

#include "itkCovariantVector.h"
#include "itkImageToImageMetric.h"
#include "itkPoint.h"

namespace anima {
template <class TFixedImage, class TMovingImage>
class FastMeanSquaresImageToImageMetric
    : public itk::ImageToImageMetric<TFixedImage, TMovingImage> {
public:
  /** Standard class typedefs. */
  using Self = FastMeanSquaresImageToImageMetric;
  using Superclass = itk::ImageToImageMetric<TFixedImage, TMovingImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(FastMeanSquaresImageToImageMetric, itk::ImageToImageMetric);

  /** Types transferred from the base class */
  using RealType = typename Superclass::RealType;
  using TransformType = typename Superclass::TransformType;
  using TransformPointer = typename Superclass::TransformPointer;
  using TransformParametersType = typename Superclass::TransformParametersType;
  using TransformJacobianType = typename Superclass::TransformJacobianType;
  using GradientPixelType = typename Superclass::GradientPixelType;
  using OutputPointType = typename Superclass::OutputPointType;
  using InputPointType = typename Superclass::InputPointType;
  using ContinuousIndexType =
      typename itk::ContinuousIndex<double, TFixedImage::ImageDimension>;

  using MeasureType = typename Superclass::MeasureType;
  using DerivativeType = typename Superclass::DerivativeType;
  using FixedImageType = typename Superclass::FixedImageType;
  using MovingImageType = typename Superclass::MovingImageType;
  using FixedImageConstPointer = typename Superclass::FixedImageConstPointer;
  using MovingImageConstPointer = typename Superclass::MovingImageConstPointer;

  /** Get the derivatives of the match measure. */
  void GetDerivative(const TransformParametersType &parameters,
                     DerivativeType &Derivative) const ITK_OVERRIDE;

  /**  Get the value for single valued optimizers. */
  MeasureType
  GetValue(const TransformParametersType &parameters) const ITK_OVERRIDE;

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative(const TransformParametersType &parameters,
                             MeasureType &Value,
                             DerivativeType &Derivative) const ITK_OVERRIDE;

  itkSetMacro(ScaleIntensities, bool);
  itkSetMacro(DefaultBackgroundValue, double);

  void PreComputeFixedValues();

protected:
  FastMeanSquaresImageToImageMetric();
  virtual ~FastMeanSquaresImageToImageMetric() {}

private:
  FastMeanSquaresImageToImageMetric(const Self &); // purposely not implemented
  void operator=(const Self &);                    // purposely not implemented

  bool m_ScaleIntensities;
  double m_DefaultBackgroundValue;

  std::vector<InputPointType> m_FixedImagePoints;
  std::vector<RealType> m_FixedImageValues;
};

} // end namespace anima

#include "animaFastMeanSquaresImageToImageMetric.hxx"
