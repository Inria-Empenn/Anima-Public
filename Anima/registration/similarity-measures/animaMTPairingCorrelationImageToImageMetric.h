#pragma once

#include <animaBaseOrientedModelImageToImageMetric.h>
#include <animaBaseTensorTools.h>
#include <animaMCMImage.h>
#include <animaMultiCompartmentModel.h>

namespace anima {

/**
 * @brief Multi-tensor correlation similarity measure as defined by Taquet et
 * al, based on pairing of the individual compartments
 *
 * M. Taquet et al. "A Mathematical Framework for the Registration and Analysis
 * of Multi-Fascicle Models for Population Studies of the Brain Microstructure".
 * IEEE TMI 2014.
 */
template <class TFixedImagePixelType, class TMovingImagePixelType,
          unsigned int ImageDimension>
class MTPairingCorrelationImageToImageMetric
    : public BaseOrientedModelImageToImageMetric<
          anima::MCMImage<TFixedImagePixelType, ImageDimension>,
          anima::MCMImage<TMovingImagePixelType, ImageDimension>> {
public:
  /** Standard class typedefs. */
  using TFixedImage = anima::MCMImage<TFixedImagePixelType, ImageDimension>;
  using TMovingImage = anima::MCMImage<TMovingImagePixelType, ImageDimension>;

  using Self = MTPairingCorrelationImageToImageMetric;
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
  itkTypeMacro(MTPairingCorrelationImageToImageMetric,
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

  void PreComputeFixedValues();

protected:
  MTPairingCorrelationImageToImageMetric();
  virtual ~MTPairingCorrelationImageToImageMetric() {}

  bool CheckTensorCompatibility() const;
  double ComputeMapping(
      const std::vector<std::vector<double>> &refImageCompartmentWeights,
      const std::vector<std::vector<PixelType>> &refImageLogTensors,
      const std::vector<std::vector<double>> &movingImageCompartmentWeights,
      const std::vector<std::vector<PixelType>> &movingImageLogTensors) const;

  bool isZero(PixelType &vector) const;

private:
  MTPairingCorrelationImageToImageMetric(
      const Self &);            // purposely not implemented
  void operator=(const Self &); // purposely not implemented

  MCModelPointer m_ZeroDiffusionModel;

  std::vector<InputPointType> m_FixedImagePoints;
  std::vector<std::vector<double>> m_FixedImageCompartmentWeights;
  std::vector<std::vector<PixelType>> m_FixedImageLogTensors;
  unsigned int m_NumberOfFixedCompartments;

  LECalculatorPointer m_leCalculator;
};

} // end namespace anima

#include "animaMTPairingCorrelationImageToImageMetric.hxx"
