#pragma once
#include <animaRefactoredBaseBMRegistrationMethod.h>

namespace anima
{

template <typename TInputImageType>
class RefactoredKissingSymmetricBMRegistrationMethod : public anima::RefactoredBaseBMRegistrationMethod <TInputImageType>
{
public:
    /** Standard class typedefs. */
    typedef RefactoredKissingSymmetricBMRegistrationMethod Self;
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
    itkTypeMacro(RefactoredKissingSymmetricBMRegistrationMethod, RefactoredBaseBMRegistrationMethod);

    itkNewMacro(Self);

protected:
    RefactoredKissingSymmetricBMRegistrationMethod() {}
    virtual ~RefactoredKissingSymmetricBMRegistrationMethod() {}

    virtual void PerformOneIteration(InputImageType *refImage, InputImageType *movingImage, TransformPointer &addOn);

private:
    RefactoredKissingSymmetricBMRegistrationMethod(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
};

} // end namespace anima

#include "animaRefactoredKissingSymmetricBMRegistrationMethod.hxx"
