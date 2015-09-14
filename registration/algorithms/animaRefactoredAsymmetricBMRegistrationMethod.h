#pragma once
#include <animaRefactoredBaseBMRegistrationMethod.h>

namespace anima
{

template <typename TInputImageType>
class RefactoredAsymmetricBMRegistrationMethod : public anima::RefactoredBaseBMRegistrationMethod <TInputImageType>
{
public:
    /** Standard class typedefs. */
    typedef RefactoredAsymmetricBMRegistrationMethod Self;
    typedef RefactoredBaseBMRegistrationMethod <TInputImageType> Superclass;
    typedef itk::SmartPointer <Self> Pointer;
    typedef itk::SmartPointer <const Self> ConstPointer;

    typedef typename Superclass::InputImageType InputImageType;
    typedef typename Superclass::TransformPointer TransformPointer;

    /** Run-time type information (and related methods). */
    itkTypeMacro(RefactoredAsymmetricBMRegistrationMethod, RefactoredBaseBMRegistrationMethod);

    itkNewMacro(Self);

protected:
    RefactoredAsymmetricBMRegistrationMethod() {}
    virtual ~RefactoredAsymmetricBMRegistrationMethod() {}

    virtual void PerformOneIteration(InputImageType *refImage, InputImageType *movingImage, TransformPointer &addOn);

private:
    RefactoredAsymmetricBMRegistrationMethod(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
};

} // end namespace anima

#include "animaRefactoredAsymmetricBMRegistrationMethod.hxx"
