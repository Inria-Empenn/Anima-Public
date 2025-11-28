#pragma once

#include <animaBaseOrientedModelImageToImageMetric.h>
#include <animaBaseTensorTools.h>
#include <animaMCMImage.h>
#include <animaMultiCompartmentModel.h>

namespace anima {

template <class TFixedImagePixelType, class TMovingImagePixelType,
          unsigned int ImageDimension>
class MCMPairingMeanSquaresImageToImageMetric
    : public BaseOrientedModelImageToImageMetric<
          anima::MCMImage<TFixedImagePixelType, ImageDimension>,
          anima::MCMImage<TMovingImagePixelType, ImageDimension>> {
public:
  /** Standard class typedefs. */
  using TFixedImage = anima::MCMImage<TFixedImagePixelType, ImageDimension>;
  using TMovingImage = anima::MCMImage<TMovingImagePixelType, ImageDimension>;

  using Self = MCMPairingMeanSquaresImageToImageMetric;
  using Superclass =
      BaseOrientedModelImageToImageMetric<TFixedImage, TMovingImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  using MCModelType = anima::MultiCompartmentModel;
  using MCModelPointer = typename MCModelType::Pointer;
  using GradientType = typename MCModelType::Vector3DType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MCMPairingMeanSquaresImageToImageMetric,
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

  using LECalculatorType = anima::LogEuclideanTensorCalculator<double>;
  using LECalculatorPointer = typename LECalculatorType::Pointer;

  /**  Get the value for single valued optimizers. */
  MeasureType
  GetValue(const TransformParametersType &parameters) const ITK_OVERRIDE;

  itkSetMacro(OneToOneMapping, bool);

  void PreComputeFixedValues();

protected:
  MCMPairingMeanSquaresImageToImageMetric();
  virtual ~MCMPairingMeanSquaresImageToImageMetric() {}

  bool CheckTensorCompatibility() const;
  double ComputeTensorBasedMetricPart(unsigned int index,
                                      const MCModelPointer &movingValue) const;
  double
  ComputeNonTensorBasedMetricPart(unsigned int index,
                                  const MCModelPointer &movingValue) const;

  bool isZero(PixelType &vector) const;

private:
  MCMPairingMeanSquaresImageToImageMetric(
      const Self &);            // purposely not implemented
  void operator=(const Self &); // purposely not implemented

  MCModelPointer m_ZeroDiffusionModel;
  PixelType m_ZeroDiffusionVector;

  std::vector<InputPointType> m_FixedImagePoints;
  std::vector<MCModelPointer> m_FixedImageValues;

  bool m_OneToOneMapping;

  LECalculatorPointer m_leCalculator;
};

} // end namespace anima

#include "animaMCMPairingMeanSquaresImageToImageMetric.hxx"
