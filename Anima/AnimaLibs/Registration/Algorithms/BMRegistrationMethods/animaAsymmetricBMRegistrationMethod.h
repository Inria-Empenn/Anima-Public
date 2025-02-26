#pragma once
#include <animaBaseBMRegistrationMethod.h>

namespace anima
{

template <typename TInputImageType>
class AsymmetricBMRegistrationMethod : public anima::BaseBMRegistrationMethod <TInputImageType>
{
public:
    /** Standard class typedefs. */
    typedef AsymmetricBMRegistrationMethod Self;
    typedef BaseBMRegistrationMethod <TInputImageType> Superclass;
    typedef itk::SmartPointer <Self> Pointer;
    typedef itk::SmartPointer <const Self> ConstPointer;

    typedef typename Superclass::InputImageType InputImageType;
    typedef typename Superclass::TransformPointer TransformPointer;
    typedef typename Superclass::AgregatorType AgregatorType;

    /** Run-time type information (and related methods). */
    itkTypeMacro(AsymmetricBMRegistrationMethod, BaseBMRegistrationMethod)

    itkNewMacro(Self)

protected:
    AsymmetricBMRegistrationMethod() {}
    virtual ~AsymmetricBMRegistrationMethod() {}

    virtual void PerformOneIteration(InputImageType *refImage, InputImageType *movingImage, TransformPointer &addOn) ITK_OVERRIDE;

private:
    AsymmetricBMRegistrationMethod(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
};

} // end namespace anima

#include "animaAsymmetricBMRegistrationMethod.hxx"
