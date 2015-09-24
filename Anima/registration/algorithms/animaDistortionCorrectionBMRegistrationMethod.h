#pragma once
#include <animaBaseBMRegistrationMethod.h>

namespace anima
{

template <typename TInputImageType>
class DistortionCorrectionBMRegistrationMethod : public anima::BaseBMRegistrationMethod <TInputImageType>
{
public:
    /** Standard class typedefs. */
    typedef DistortionCorrectionBMRegistrationMethod Self;
    typedef BaseBMRegistrationMethod <TInputImageType> Superclass;
    typedef itk::SmartPointer <Self> Pointer;
    typedef itk::SmartPointer <const Self> ConstPointer;

    typedef typename Superclass::ImageScalarType ImageScalarType;
    typedef typename Superclass::InputImageType InputImageType;
    typedef typename Superclass::InputImagePointer InputImagePointer;
    typedef typename Superclass::TransformType TransformType;
    typedef typename Superclass::TransformPointer TransformPointer;
    typedef typename Superclass::BlockMatcherType BlockMatcherType;
    typedef typename Superclass::AgregatorType AgregatorType;
    typedef typename Superclass::AgregatorScalarType AgregatorScalarType;
    typedef typename Superclass::AffineTransformType AffineTransformType;
    typedef typename Superclass::SVFTransformType SVFTransformType;
    typedef typename Superclass::DisplacementFieldTransformType DisplacementFieldTransformType;
    typedef typename Superclass::DisplacementFieldTransformPointer DisplacementFieldTransformPointer;

    /** Run-time type information (and related methods). */
    itkTypeMacro(DistortionCorrectionBMRegistrationMethod, BaseBMRegistrationMethod);

    itkNewMacro(Self);

    itkSetMacro(CurrentTransform, TransformPointer);

protected:
    DistortionCorrectionBMRegistrationMethod() {}
    virtual ~DistortionCorrectionBMRegistrationMethod() {}

    virtual void SetupTransform(TransformPointer &optimizedTransform);
    virtual void PerformOneIteration(InputImageType *refImage, InputImageType *movingImage, TransformPointer &addOn);
    virtual void ResampleImages(TransformType *currentTransform, InputImagePointer &refImage, InputImagePointer &movingImage);
    virtual bool ComposeAddOnWithTransform(TransformPointer &computedTransform, TransformType *addOn);

private:
    DistortionCorrectionBMRegistrationMethod(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    TransformPointer m_CurrentTransform;
};

} // end namespace anima

#include "animaDistortionCorrectionBMRegistrationMethod.hxx"
