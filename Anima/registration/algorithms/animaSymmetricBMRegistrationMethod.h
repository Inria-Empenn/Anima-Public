#pragma once
#include <animaBaseBMRegistrationMethod.h>

namespace anima {

template <typename TInputImageType>
class SymmetricBMRegistrationMethod
    : public anima::BaseBMRegistrationMethod<TInputImageType> {
public:
  /** Standard class typedefs. */
  using Self = SymmetricBMRegistrationMethod;
  using Superclass = BaseBMRegistrationMethod<TInputImageType>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  using InputImageType = typename Superclass::InputImageType;
  using TransformPointer = typename Superclass::TransformPointer;
  using BlockMatcherType = typename Superclass::BlockMatcherType;
  using AgregatorType = typename Superclass::AgregatorType;
  using AffineTransformType = typename Superclass::AffineTransformType;
  using SVFTransformType = typename Superclass::SVFTransformType;

  /** Run-time type information (and related methods). */
  itkTypeMacro(SymmetricBMRegistrationMethod, BaseBMRegistrationMethod);

  itkNewMacro(Self);

  void SetReverseBlockMatcher(BlockMatcherType *matcher) {
    m_ReverseBlockMatcher = matcher;
  }

protected:
  SymmetricBMRegistrationMethod() {}
  virtual ~SymmetricBMRegistrationMethod() {}

  virtual void PerformOneIteration(InputImageType *refImage,
                                   InputImageType *movingImage,
                                   TransformPointer &addOn) ITK_OVERRIDE;

private:
  SymmetricBMRegistrationMethod(const Self &); // purposely not implemented
  void operator=(const Self &);                // purposely not implemented

  BlockMatcherType *m_ReverseBlockMatcher;
};

} // end namespace anima

#include "animaSymmetricBMRegistrationMethod.hxx"
