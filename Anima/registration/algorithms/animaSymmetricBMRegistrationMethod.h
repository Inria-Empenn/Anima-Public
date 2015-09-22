#pragma once
#include <animaBaseBMRegistrationMethod.h>

namespace anima
{

template <typename TInputImageType>
class SymmetricBMRegistrationMethod : public anima::BaseBMRegistrationMethod <TInputImageType>
{
public:
    /** Standard class typedefs. */
    typedef SymmetricBMRegistrationMethod Self;
    typedef BaseBMRegistrationMethod <TInputImageType> Superclass;
    typedef itk::SmartPointer <Self> Pointer;
    typedef itk::SmartPointer <const Self> ConstPointer;

    typedef typename Superclass::InputImageType InputImageType;
    typedef typename Superclass::TransformPointer TransformPointer;
    typedef typename Superclass::BlockMatcherType BlockMatcherType;
    typedef typename Superclass::AgregatorType AgregatorType;
    typedef typename Superclass::AffineTransformType AffineTransformType;
    typedef typename Superclass::SVFTransformType SVFTransformType;

    /** Run-time type information (and related methods). */
    itkTypeMacro(SymmetricBMRegistrationMethod, BaseBMRegistrationMethod);

    itkNewMacro(Self);

    void SetReverseBlockMatcher(BlockMatcherType *matcher) {m_ReverseBlockMatcher = matcher;}

protected:
    SymmetricBMRegistrationMethod() {}
    virtual ~SymmetricBMRegistrationMethod() {}

    virtual void PerformOneIteration(InputImageType *refImage, InputImageType *movingImage, TransformPointer &addOn);

private:
    SymmetricBMRegistrationMethod(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    BlockMatcherType *m_ReverseBlockMatcher;
};

} // end namespace anima

#include "animaSymmetricBMRegistrationMethod.hxx"
