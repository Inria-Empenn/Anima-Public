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
    typedef typename Superclass::TransformType TransformType;
    typedef typename Superclass::TransformPointer TransformPointer;
    typedef typename Superclass::BlockMatcherType BlockMatcherType;
    typedef typename Superclass::AgregatorType AgregatorType;
    typedef typename Superclass::AffineTransformType AffineTransformType;
    typedef typename Superclass::AffineTransformPointer AffineTransformPointer;
    typedef typename Superclass::DisplacementFieldTransformType DisplacementFieldTransformType;
    typedef typename Superclass::DisplacementFieldTransformPointer DisplacementFieldTransformPointer;
    typedef typename Superclass::SVFTransformType SVFTransformType;

    /** Run-time type information (and related methods). */
    itkTypeMacro(KissingSymmetricBMRegistrationMethod, BaseBMRegistrationMethod)

    itkNewMacro(Self)

    itkSetMacro(ReferenceBackgroundValue, double)
    itkSetMacro(FloatingBackgroundValue, double)
    itkSetMacro(RegistrationPointLocation, double)

protected:
    KissingSymmetricBMRegistrationMethod()
    {
        m_ReferenceBackgroundValue = 0;
        m_FloatingBackgroundValue = 0;
        m_RegistrationPointLocation = 0.5;
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

    virtual TransformPointer GetBackwardTransformForResampling(TransformType *transform) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(KissingSymmetricBMRegistrationMethod);

    double m_ReferenceBackgroundValue;
    double m_FloatingBackgroundValue;

    //! Registers the two images towards a point located in the range [0, 1]: 0 denotes on ref, 1: on moving, anything else lies on the path
    double m_RegistrationPointLocation;
};

} // end namespace anima

#include "animaKissingSymmetricBMRegistrationMethod.hxx"
