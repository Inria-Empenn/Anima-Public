#pragma once
#include <animaBaseBMRegistrationMethod.h>

namespace anima
{

template <typename TInputImageType>
class KissingSymmetricBMRegistrationMethod : public anima::BaseBMRegistrationMethod <TInputImageType>
{
public:
    /** Standard class typedefs. */
    typedef KissingSymmetricBMRegistrationMethod Self;
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
    itkTypeMacro(KissingSymmetricBMRegistrationMethod, BaseBMRegistrationMethod)

    itkNewMacro(Self)

protected:
    KissingSymmetricBMRegistrationMethod() {}
    virtual ~KissingSymmetricBMRegistrationMethod() {}

    virtual void PerformOneIteration(InputImageType *refImage, InputImageType *movingImage, TransformPointer &addOn) ITK_OVERRIDE;

private:
    KissingSymmetricBMRegistrationMethod(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
};

} // end namespace anima

#include "animaKissingSymmetricBMRegistrationMethod.hxx"
