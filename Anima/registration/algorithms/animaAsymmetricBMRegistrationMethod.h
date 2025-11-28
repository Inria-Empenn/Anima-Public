#pragma once
#include <animaBaseBMRegistrationMethod.h>

namespace anima {

template <typename TInputImageType>
class AsymmetricBMRegistrationMethod
    : public anima::BaseBMRegistrationMethod<TInputImageType> {
public:
  /** Standard class typedefs. */
  using Self = AsymmetricBMRegistrationMethod;
  using Superclass = BaseBMRegistrationMethod<TInputImageType>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  using InputImageType = typename Superclass::InputImageType;
  using TransformPointer = typename Superclass::TransformPointer;
  using AgregatorType = typename Superclass::AgregatorType;

  /** Run-time type information (and related methods). */
  itkTypeMacro(AsymmetricBMRegistrationMethod, BaseBMRegistrationMethod);

  itkNewMacro(Self);

protected:
  AsymmetricBMRegistrationMethod() {}
  virtual ~AsymmetricBMRegistrationMethod() {}

  virtual void PerformOneIteration(InputImageType *refImage,
                                   InputImageType *movingImage,
                                   TransformPointer &addOn) ITK_OVERRIDE;

private:
  AsymmetricBMRegistrationMethod(const Self &); // purposely not implemented
  void operator=(const Self &);                 // purposely not implemented
};

} // end namespace anima

#include "animaAsymmetricBMRegistrationMethod.hxx"
