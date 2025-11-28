#pragma once

#include <animaBaseOrientedModelImageToImageMetric.h>
#include <animaMCMImage.h>
#include <animaMultiCompartmentModel.h>

namespace anima {

template <class TFixedImagePixelType, class TMovingImagePixelType,
          unsigned int ImageDimension>
class MCMCorrelationImageToImageMetric
    : public BaseOrientedModelImageToImageMetric<
          anima::MCMImage<TFixedImagePixelType, ImageDimension>,
          anima::MCMImage<TMovingImagePixelType, ImageDimension>> {
public:
  /** Standard class typedefs. */
  using TFixedImage = anima::MCMImage<TFixedImagePixelType, ImageDimension>;
  using TMovingImage = anima::MCMImage<TMovingImagePixelType, ImageDimension>;

  using Self = MCMCorrelationImageToImageMetric;
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
  itkTypeMacro(MCMCorrelationImageToImageMetric,
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

  void SetSmallDelta(double val);
  void SetBigDelta(double val);
  void SetGradientStrengths(std::vector<double> &val);
  void SetGradientDirections(std::vector<GradientType> &val);

  itkSetMacro(ForceApproximation, bool);
  itkSetMacro(LowerBoundGaussianSigma, double);
  itkSetMacro(UpperBoundGaussianSigma, double);

  void PreComputeFixedValues();

protected:
  MCMCorrelationImageToImageMetric();
  virtual ~MCMCorrelationImageToImageMetric() {}

  //! Compute base integration weights for non tensor integration
  void UpdateSphereWeights();

  bool CheckTensorCompatibility() const;
  double ComputeTensorBasedMetric(
      const std::vector<MCModelPointer> &movingValues) const;
  double ComputeNonTensorBasedMetric(
      const std::vector<MCModelPointer> &movingValues) const;

  bool isZero(PixelType &vector) const;

private:
  MCMCorrelationImageToImageMetric(const Self &); // purposely not implemented
  void operator=(const Self &);                   // purposely not implemented

  std::vector<InputPointType> m_FixedImagePoints;
  std::vector<MCModelPointer> m_FixedImageValues;

  // Optional parameters for the case when compartments are not tensor
  // compatible
  std::vector<double> m_GradientStrengths;
  double m_SmallDelta, m_BigDelta;
  std::vector<GradientType> m_GradientDirections;

  // Parameters for numerical integration on non tensor compatible models
  std::vector<double> m_SphereWeights;
  std::vector<unsigned int> m_BValWeightsIndexes;

  MCModelPointer m_ZeroDiffusionModel;
  PixelType m_ZeroDiffusionVector;

  bool m_ForceApproximation;

  // Lower and upper bounds of Gaussian sigma for smoothing
  double m_LowerBoundGaussianSigma;
  double m_UpperBoundGaussianSigma;
};

} // end namespace anima

#include "animaMCMCorrelationImageToImageMetric.hxx"
