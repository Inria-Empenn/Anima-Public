#pragma once

#include <animaBaseOrientedModelImageToImageMetric.h>
#include <animaMCMImage.h>
#include <animaMCML2DistanceComputer.h>
#include <animaMultiCompartmentModel.h>

namespace anima {

template <class TFixedImagePixelType, class TMovingImagePixelType,
          unsigned int ImageDimension>
class MCMMeanSquaresImageToImageMetric
    : public BaseOrientedModelImageToImageMetric<
          anima::MCMImage<TFixedImagePixelType, ImageDimension>,
          anima::MCMImage<TMovingImagePixelType, ImageDimension>> {
public:
  /** Standard class typedefs. */
  using TFixedImage = anima::MCMImage<TFixedImagePixelType, ImageDimension>;
  using TMovingImage = anima::MCMImage<TMovingImagePixelType, ImageDimension>;

  using Self = MCMMeanSquaresImageToImageMetric;
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
  itkTypeMacro(MCMMeanSquaresImageToImageMetric,
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

  using MCML2DistanceComputerPointer = anima::MCML2DistanceComputer::Pointer;

  /**  Get the value for single valued optimizers. */
  MeasureType
  GetValue(const TransformParametersType &parameters) const ITK_OVERRIDE;

  void SetSmallDelta(double val) { m_L2DistanceComputer->SetSmallDelta(val); }
  void SetBigDelta(double val) { m_L2DistanceComputer->SetBigDelta(val); }
  void SetGradientStrengths(std::vector<double> &val) {
    m_L2DistanceComputer->SetGradientStrengths(val);
  }
  void SetGradientDirections(std::vector<GradientType> &val) {
    m_L2DistanceComputer->SetGradientDirections(val);
  }

  void SetForceApproximation(bool val) {
    m_L2DistanceComputer->SetForceApproximation(val);
  }
  void SetLowPassGaussianSigma(double val) {
    m_L2DistanceComputer->SetLowPassGaussianSigma(val);
  }

  void PreComputeFixedValues();

protected:
  MCMMeanSquaresImageToImageMetric();
  virtual ~MCMMeanSquaresImageToImageMetric() {}

  bool isZero(PixelType &vector) const;

private:
  MCMMeanSquaresImageToImageMetric(const Self &); // purposely not implemented
  void operator=(const Self &);                   // purposely not implemented

  std::vector<InputPointType> m_FixedImagePoints;
  std::vector<MCModelPointer> m_FixedImageValues;

  MCML2DistanceComputerPointer m_L2DistanceComputer;

  MCModelPointer m_ZeroDiffusionModel;
  PixelType m_ZeroDiffusionVector;
};

} // end namespace anima

#include "animaMCMMeanSquaresImageToImageMetric.hxx"
