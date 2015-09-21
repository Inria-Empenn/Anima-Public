#pragma once
#include <animaRefactoredBaseBMRegistrationMethod.h>

namespace anima
{

template <typename TInputImageType>
class RefactoredSymmetricBMRegistrationMethod : public anima::RefactoredBaseBMRegistrationMethod <TInputImageType>
{
public:
    /** Standard class typedefs. */
    typedef RefactoredSymmetricBMRegistrationMethod Self;
    typedef RefactoredBaseBMRegistrationMethod <TInputImageType> Superclass;
    typedef itk::SmartPointer <Self> Pointer;
    typedef itk::SmartPointer <const Self> ConstPointer;

    typedef typename Superclass::InputImageType InputImageType;
    typedef typename Superclass::TransformPointer TransformPointer;
    typedef typename Superclass::BlockMatcherType BlockMatcherType;
    typedef typename Superclass::AgregatorType AgregatorType;
    typedef typename Superclass::AffineTransformType AffineTransformType;
    typedef typename Superclass::SVFTransformType SVFTransformType;

    /** Run-time type information (and related methods). */
    itkTypeMacro(RefactoredSymmetricBMRegistrationMethod, RefactoredBaseBMRegistrationMethod);

    itkNewMacro(Self);

    void SetReverseBlockMatcher(BlockMatcherType *matcher) {m_ReverseBlockMatcher = matcher;}

protected:
    RefactoredSymmetricBMRegistrationMethod() {}
    virtual ~RefactoredSymmetricBMRegistrationMethod() {}

    virtual void PerformOneIteration(InputImageType *refImage, InputImageType *movingImage, TransformPointer &addOn);

private:
    RefactoredSymmetricBMRegistrationMethod(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    BlockMatcherType *m_ReverseBlockMatcher;
};

} // end namespace anima

#include "animaRefactoredSymmetricBMRegistrationMethod.hxx"
