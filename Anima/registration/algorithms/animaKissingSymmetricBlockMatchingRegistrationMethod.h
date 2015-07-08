#pragma once

#include <animaBlockMatchingBaseRegistrationMethod.h>
#include <animaBlockMatchInitializer.h>

namespace anima
{

template <typename TInputImage>
class KissingSymmetricBlockMatchingRegistrationMethod : public BlockMatchingBaseRegistrationMethod <TInputImage>
{
public:
    /** Standard class typedefs. */
    typedef KissingSymmetricBlockMatchingRegistrationMethod Self;
    typedef BlockMatchingBaseRegistrationMethod <TInputImage> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    typedef typename Superclass::InputImageType InputImageType;
    typedef typename Superclass::TransformType TransformType;
    typedef typename Superclass::TransformPointer TransformPointer;
    typedef typename Superclass::ThreadedMatchData ThreadedMatchData;
    typedef typename Superclass::ImageIndexType ImageIndexType;
    typedef typename Superclass::ImageRegionType ImageRegionType;
    typedef typename Superclass::AffineTransformType AffineTransformType;
    typedef typename Superclass::SVFTransformType SVFTransformType;
    typedef typename Superclass::AgregatorType AgregatorType;

    typedef typename TInputImage::PixelType InputPixelType;
    typedef typename anima::BlockMatchingInitializer<InputPixelType,TInputImage::ImageDimension> InitializerType;
    typedef typename InitializerType::Pointer InitializerPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(KissingSymmetricBlockMatchingRegistrationMethod, BlockMatchingBaseRegistrationMethod);

    itkSetMacro (BlockPercentageKept, double);
    itkSetMacro (BlockSize, unsigned int);
    itkSetMacro (BlockSpacing, unsigned int);

    itkSetMacro (UseTransformationDam, bool);
    itkSetMacro (DamDistance, double);

protected:
    KissingSymmetricBlockMatchingRegistrationMethod();
    virtual ~KissingSymmetricBlockMatchingRegistrationMethod() {}

    virtual void PerformOneIteration(InputImageType *refImage, InputImageType *movingImage, TransformPointer &addOn);

    void InitializeBlocksOnImage(InitializerPointer &initPtr, InputImageType *image);

private:
    KissingSymmetricBlockMatchingRegistrationMethod(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    // Block generation parameters
    double m_BlockPercentageKept;
    unsigned int m_BlockSize;
    unsigned int m_BlockSpacing;
    bool m_UseTransformationDam;
    double m_DamDistance;
};

} // end namespace anima

#include "animaKissingSymmetricBlockMatchingRegistrationMethod.hxx"
