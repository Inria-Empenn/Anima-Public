#pragma once

#include <itkImageBase.h>
#include <itkInterpolateImageFunction.h>
#include <itkMacro.h>
#include <itkSingleValuedCostFunction.h>
#include <itkSpatialObject.h>
#include <itkTransform.h>

namespace anima {
template <class TFixedImage, class TMovingImage>
class BaseOrientedModelImageToImageMetric
    : public itk::SingleValuedCostFunction {
public:
  /** Standard class typedefs. */
  using Self = BaseOrientedModelImageToImageMetric;
  using Superclass = itk::SingleValuedCostFunction;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Type used for representing point components  */
  using CoordinateRepresentationType = typename Superclass::ParametersValueType;

  /** Run-time type information (and related methods). */
  itkTypeMacro(BaseOrientedModelImageToImageMetric, SingleValuedCostFunction);

  /**  Type of the moving Image. */
  using MovingImageType = TMovingImage;
  using MovingImagePixelType = typename TMovingImage::PixelType;
  using MovingImageConstPointer = typename MovingImageType::ConstPointer;

  /**  Type of the fixed Image. */
  using FixedImageType = TFixedImage;
  using FixedImageConstPointer = typename FixedImageType::ConstPointer;
  using FixedImageRegionType = typename FixedImageType::RegionType;

  /** Constants for the image dimensions */
  itkStaticConstMacro(MovingImageDimension, unsigned int,
                      TMovingImage::ImageDimension);
  itkStaticConstMacro(FixedImageDimension, unsigned int,
                      TFixedImage::ImageDimension);

  /**  Type of the Transform Base class */
  using TransformType =
      itk::Transform<CoordinateRepresentationType,
                     itkGetStaticConstMacro(MovingImageDimension),
                     itkGetStaticConstMacro(FixedImageDimension)>;

  using TransformPointer = typename TransformType::Pointer;
  using InputPointType = typename TransformType::InputPointType;
  using OutputPointType = typename TransformType::OutputPointType;
  using TransformParametersType = typename TransformType::ParametersType;
  using TransformJacobianType = typename TransformType::JacobianType;

  /**  Type of the Interpolator Base class */
  using InterpolatorType =
      itk::InterpolateImageFunction<MovingImageType,
                                    CoordinateRepresentationType>;

  using InterpolatorPointer = typename InterpolatorType::Pointer;

  /**  Type for the mask of the fixed image. Only pixels that are "inside"
   this mask will be considered for the computation of the metric */
  using FixedImageMaskType =
      itk::SpatialObject<itkGetStaticConstMacro(FixedImageDimension)>;
  using FixedImageMaskPointer = typename FixedImageMaskType::Pointer;
  using FixedImageMaskConstPointer = typename FixedImageMaskType::ConstPointer;

  /**  Type for the mask of the moving image. Only pixels that are "inside"
   this mask will be considered for the computation of the metric */
  using MovingImageMaskType =
      itk::SpatialObject<itkGetStaticConstMacro(MovingImageDimension)>;
  using MovingImageMaskPointer = typename MovingImageMaskType::Pointer;
  using MovingImageMaskConstPointer =
      typename MovingImageMaskType::ConstPointer;

  using RealType = double;

  /**  Type of the measure. */
  using MeasureType = typename Superclass::MeasureType;

  /**  Type of the derivative. */
  using DerivativeType = typename Superclass::DerivativeType;

  /**  Type of the parameters. */
  using ParametersType = typename Superclass::ParametersType;

  /** Connect the Fixed Image.  */
  itkSetConstObjectMacro(FixedImage, FixedImageType);

  /** Get the Fixed Image. */
  itkGetConstObjectMacro(FixedImage, FixedImageType);

  /** Connect the Moving Image.  */
  itkSetConstObjectMacro(MovingImage, MovingImageType);

  /** Get the Moving Image. */
  itkGetConstObjectMacro(MovingImage, MovingImageType);

  /** Connect the Transform. */
  itkSetObjectMacro(Transform, TransformType);

  /** Get a pointer to the Transform.  */
  itkGetConstObjectMacro(Transform, TransformType);

  /** Connect the Interpolator. */
  itkSetObjectMacro(Interpolator, InterpolatorType);

  /** Get a pointer to the Interpolator.  */
  itkGetConstObjectMacro(Interpolator, InterpolatorType);

  /** Get the number of pixels considered in the computation. */
  itkGetConstReferenceMacro(NumberOfPixelsCounted, unsigned long);

  /** Set the region over which the metric will be computed */
  itkSetMacro(FixedImageRegion, FixedImageRegionType);

  /** Get the region over which the metric will be computed */
  itkGetConstReferenceMacro(FixedImageRegion, FixedImageRegionType);

  /** Set/Get the moving image mask. */
  itkSetObjectMacro(MovingImageMask, MovingImageMaskType);
  itkSetConstObjectMacro(MovingImageMask, MovingImageMaskType);
  itkGetConstObjectMacro(MovingImageMask, MovingImageMaskType);

  /** Set/Get the fixed image mask. */
  itkSetObjectMacro(FixedImageMask, FixedImageMaskType);
  itkSetConstObjectMacro(FixedImageMask, FixedImageMaskType);
  itkGetConstObjectMacro(FixedImageMask, FixedImageMaskType);

  /** Set the parameters defining the Transform. */
  void SetTransformParameters(const ParametersType &parameters) const;

  /** Return the number of parameters required by the Transform */
  unsigned int GetNumberOfParameters(void) const ITK_OVERRIDE {
    return m_Transform->GetNumberOfParameters();
  }

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  virtual void Initialize();

  //! Should not be used
  void GetDerivative(const ParametersType &parameters,
                     DerivativeType &derivative) const ITK_OVERRIDE {
    itkExceptionMacro("No derivatives implemented for tensor metrics...");
  }

  void GetValueAndDerivative(const TransformParametersType &parameters,
                             MeasureType &value,
                             DerivativeType &derivative) const ITK_OVERRIDE {
    itkExceptionMacro("No derivatives implemented for tensor metrics...");
  }

  enum ModelReorientationType { NONE = 0, FINITE_STRAIN, PPD };

  itkSetMacro(ModelRotation, ModelReorientationType);
  itkGetConstMacro(ModelRotation, ModelReorientationType);

protected:
  BaseOrientedModelImageToImageMetric();
  virtual ~BaseOrientedModelImageToImageMetric();
  void PrintSelf(std::ostream &os, itk::Indent indent) const ITK_OVERRIDE;

  mutable unsigned long m_NumberOfPixelsCounted;

  FixedImageConstPointer m_FixedImage;
  MovingImageConstPointer m_MovingImage;

  mutable TransformPointer m_Transform;
  mutable vnl_matrix<double> m_OrientationMatrix;
  InterpolatorPointer m_Interpolator;

  FixedImageMaskConstPointer m_FixedImageMask;
  MovingImageMaskConstPointer m_MovingImageMask;

private:
  BaseOrientedModelImageToImageMetric(const Self &); // purposely not
                                                     // implemented
  void operator=(const Self &); // purposely not implemented

  FixedImageRegionType m_FixedImageRegion;

  ModelReorientationType m_ModelRotation;
};

} // end of namespace anima

#include "animaBaseOrientedModelImageToImageMetric.hxx"
