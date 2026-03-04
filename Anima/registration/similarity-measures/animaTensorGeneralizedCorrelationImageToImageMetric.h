#pragma once

#include <animaBaseOrientedModelImageToImageMetric.h>
#include <animaBaseTensorTools.h>

#include <itkCovariantVector.h>
#include <itkPoint.h>
#include <itkVectorImage.h>

namespace anima {
template <class TFixedImagePixelType, class TMovingImagePixelType,
          unsigned int ImageDimension>
class TensorGeneralizedCorrelationImageToImageMetric
    : public anima::BaseOrientedModelImageToImageMetric<
          itk::VectorImage<TFixedImagePixelType, ImageDimension>,
          itk::VectorImage<TMovingImagePixelType, ImageDimension>> {
public:
  /** Standard class typedefs. */
  using TFixedImage = itk::VectorImage<TFixedImagePixelType, ImageDimension>;
  using TMovingImage = itk::VectorImage<TMovingImagePixelType, ImageDimension>;

  using Self = TensorGeneralizedCorrelationImageToImageMetric;
  using Superclass =
      anima::BaseOrientedModelImageToImageMetric<TFixedImage, TMovingImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(TensorGeneralizedCorrelationImageToImageMetric,
               BaseOrientedModelImageToImageMetric);

  /** Types transferred from the base class */
  using PixelType = typename TFixedImage::PixelType;
  using CovarianceType = vnl_matrix<double>;

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

  using LECalculatorType = anima::LogEuclideanTensorCalculator<double>;
  using LECalculatorPointer = LECalculatorType::Pointer;

  /**  Get the value for single valued optimizers. */
  MeasureType
  GetValue(const TransformParametersType &parameters) const ITK_OVERRIDE;

  void PreComputeFixedValues();

  itkSetMacro(OrientationPenalty, bool);
  itkSetMacro(VarianceThreshold, double);

protected:
  TensorGeneralizedCorrelationImageToImageMetric();
  virtual ~TensorGeneralizedCorrelationImageToImageMetric() {}
  void PrintSelf(std::ostream &os, itk::Indent indent) const ITK_OVERRIDE;

private:
  TensorGeneralizedCorrelationImageToImageMetric(
      const Self &);            // purposely not implemented
  void operator=(const Self &); // purposely not implemented

  bool m_OrientationPenalty;

  double m_VarianceThreshold;

  PixelType m_FixedMean;
  CovarianceType m_FixedHalfInvCovarianceMatrix;

  std::vector<InputPointType> m_FixedImagePoints;
  std::vector<PixelType> m_FixedImageValues;

  LECalculatorPointer m_leCalculator;
};

} // end namespace anima

#include "animaTensorGeneralizedCorrelationImageToImageMetric.hxx"
