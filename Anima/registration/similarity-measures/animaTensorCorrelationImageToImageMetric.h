#pragma once

#include <animaBaseOrientedModelImageToImageMetric.h>
#include <itkCovariantVector.h>
#include <itkPoint.h>
#include <itkVectorImage.h>

namespace anima {

/**
 * @brief Tensor correlation similarity measure as defined by Taquet et al
 *
 * M. Taquet et al. "A Mathematical Framework for the Registration and Analysis
 * of Multi-Fascicle Models for Population Studies of the Brain Microstructure".
 * IEEE TMI 2014.
 */
template <class TFixedImagePixelType, class TMovingImagePixelType,
          unsigned int ImageDimension>
class TensorCorrelationImageToImageMetric
    : public BaseOrientedModelImageToImageMetric<
          itk::VectorImage<TFixedImagePixelType, ImageDimension>,
          itk::VectorImage<TMovingImagePixelType, ImageDimension>> {
public:
  /** Standard class typedefs. */
  using TFixedImage = itk::VectorImage<TFixedImagePixelType, ImageDimension>;
  using TMovingImage = itk::VectorImage<TMovingImagePixelType, ImageDimension>;

  using Self = TensorCorrelationImageToImageMetric;
  using Superclass =
      BaseOrientedModelImageToImageMetric<TFixedImage, TMovingImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(TensorCorrelationImageToImageMetric,
               BaseOrientedModelImageToImageMetric);

  /** Types transferred from the base class */
  using PixelType = typename TFixedImage::PixelType;

  using TransformType = typename Superclass::TransformType;
  using TransformPointer = typename Superclass::TransformPointer;
  using TransformParametersType = typename Superclass::TransformParametersType;
  using OutputPointType = typename Superclass::OutputPointType;
  using InputPointType = typename Superclass::InputPointType;
  using ContinuousIndexType =
      typename itk::ContinuousIndex<double, ImageDimension>;

  using CoordinateRepresentationType =
      typename Superclass::CoordinateRepresentationType;

  using MeasureType = typename Superclass::MeasureType;
  using FixedImageType = typename Superclass::FixedImageType;
  using MovingImageType = typename Superclass::MovingImageType;
  using FixedImageConstPointer = typename Superclass::FixedImageConstPointer;
  using MovingImageConstPointer = typename Superclass::MovingImageConstPointer;

  /**  Get the value for single valued optimizers. */
  MeasureType
  GetValue(const TransformParametersType &parameters) const ITK_OVERRIDE;

  void PreComputeFixedValues();

protected:
  TensorCorrelationImageToImageMetric();
  virtual ~TensorCorrelationImageToImageMetric() {}

private:
  TensorCorrelationImageToImageMetric(const Self &); // purposely not
                                                     // implemented
  void operator=(const Self &); // purposely not implemented

  double m_FixedTProduct, m_FixedDenominator;
  double m_LogEpsilon;

  std::vector<InputPointType> m_FixedImagePoints;
  std::vector<PixelType> m_FixedImageValues;
};

} // end namespace anima

#include "animaTensorCorrelationImageToImageMetric.hxx"
