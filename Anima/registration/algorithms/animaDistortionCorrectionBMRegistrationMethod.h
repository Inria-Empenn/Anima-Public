#pragma once
#include <animaBaseBMRegistrationMethod.h>

namespace anima {

template <typename TInputImageType>
class DistortionCorrectionBMRegistrationMethod
    : public anima::BaseBMRegistrationMethod<TInputImageType> {
public:
  /** Standard class typedefs. */
  using Self = DistortionCorrectionBMRegistrationMethod;
  using Superclass = BaseBMRegistrationMethod<TInputImageType>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  using ImageScalarType = typename Superclass::ImageScalarType;
  using InputImageType = typename Superclass::InputImageType;
  using InputImagePointer = typename Superclass::InputImagePointer;
  using TransformType = typename Superclass::TransformType;
  using TransformPointer = typename Superclass::TransformPointer;
  using BlockMatcherType = typename Superclass::BlockMatcherType;
  using AgregatorType = typename Superclass::AgregatorType;
  using AgregatorScalarType = typename Superclass::AgregatorScalarType;
  using AffineTransformType = typename Superclass::AffineTransformType;
  using SVFTransformType = typename Superclass::SVFTransformType;
  using DisplacementFieldTransformType =
      typename Superclass::DisplacementFieldTransformType;
  using DisplacementFieldTransformPointer =
      typename Superclass::DisplacementFieldTransformPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(DistortionCorrectionBMRegistrationMethod,
               BaseBMRegistrationMethod);

  itkNewMacro(Self);

  itkSetMacro(CurrentTransform, TransformPointer);

protected:
  DistortionCorrectionBMRegistrationMethod() {}
  virtual ~DistortionCorrectionBMRegistrationMethod() {}

  virtual void
  SetupTransform(TransformPointer &optimizedTransform) ITK_OVERRIDE;
  virtual void PerformOneIteration(InputImageType *refImage,
                                   InputImageType *movingImage,
                                   TransformPointer &addOn) ITK_OVERRIDE;
  virtual void ResampleImages(TransformType *currentTransform,
                              InputImagePointer &refImage,
                              InputImagePointer &movingImage) ITK_OVERRIDE;
  virtual bool ComposeAddOnWithTransform(TransformPointer &computedTransform,
                                         TransformType *addOn) ITK_OVERRIDE;

private:
  DistortionCorrectionBMRegistrationMethod(
      const Self &);            // purposely not implemented
  void operator=(const Self &); // purposely not implemented

  TransformPointer m_CurrentTransform;
};

} // end namespace anima

#include "animaDistortionCorrectionBMRegistrationMethod.hxx"
