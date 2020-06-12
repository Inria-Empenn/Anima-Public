#pragma once

#include <animaBaseBMRegistrationMethod.h>
#include <animaAnatomicalBlockMatcher.h>

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

    itkSetMacro(ReferenceBackgroundValue, double)
    itkSetMacro(FloatingBackgroundValue, double)

protected:
    KissingSymmetricBMRegistrationMethod()
    {
        m_ReferenceBackgroundValue = 0;
        m_FloatingBackgroundValue = 0;
    }

    virtual ~KissingSymmetricBMRegistrationMethod() {}

    virtual void PerformOneIteration(InputImageType *refImage, InputImageType *movingImage, TransformPointer &addOn) ITK_OVERRIDE;

    template <class ScalarPixelType> void ChangeDefaultBackgroundValue(itk::Image <ScalarPixelType, InputImageType::ImageDimension> *itkNotUsed(image), double backgroundValue)
    {
        using BMType = anima::AnatomicalBlockMatcher <InputImageType>;
        BMType *tmpBM = dynamic_cast <BMType *> (this->GetBlockMatcher());
        tmpBM->SetDefaultBackgroundValue(backgroundValue);
    }

    template <class ScalarPixelType> void ChangeDefaultBackgroundValue(itk::VectorImage <ScalarPixelType, InputImageType::ImageDimension> *itkNotUsed(image), double itkNotUsed(backgroundValue))
    {
        // Vector image : do nothing
    }

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(KissingSymmetricBMRegistrationMethod);

    double m_ReferenceBackgroundValue;
    double m_FloatingBackgroundValue;
};

} // end namespace anima

#include "animaKissingSymmetricBMRegistrationMethod.hxx"
