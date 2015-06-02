#pragma once

#include <animaTensorBlockMatchingBaseRegistrationMethod.h>
#include <animaBlockMatchInitializer.h>

namespace anima
{

    template <typename TInputImage>
    class TensorKissingSymmetricBlockMatchingRegistrationMethod : public TensorBlockMatchingBaseRegistrationMethod <TInputImage>
    {
    public:
        /** Standard class typedefs. */
        typedef TensorKissingSymmetricBlockMatchingRegistrationMethod Self;
        typedef TensorBlockMatchingBaseRegistrationMethod <TInputImage> Superclass;
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

        typedef typename TInputImage::IOPixelType InputPixelType;
        typedef typename anima::BlockMatchingInitializer<InputPixelType,TInputImage::ImageDimension> InitializerType;
        typedef typename InitializerType::Pointer InitializerPointer;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro(TensorKissingSymmetricBlockMatchingRegistrationMethod, BlockMatchingBaseRegistrationMethod);

        itkSetMacro (BlockPercentageKept, double);
        itkSetMacro (BlockSize, unsigned int);
        itkSetMacro (BlockSpacing, unsigned int);

        itkSetMacro (UseTransformationDam, bool);
        itkSetMacro (DamDistance, double);

    protected:
        TensorKissingSymmetricBlockMatchingRegistrationMethod();
        virtual ~TensorKissingSymmetricBlockMatchingRegistrationMethod() {}

        virtual void PerformOneIteration(InputImageType *refImage, InputImageType *movingImage, TransformPointer &addOn);

        void InitializeBlocksOnImage(InitializerPointer &initPtr, InputImageType *image);

    private:
        TensorKissingSymmetricBlockMatchingRegistrationMethod(const Self&); //purposely not implemented
        void operator=(const Self&); //purposely not implemented

        // Block generation parameters
        double m_BlockPercentageKept;
        unsigned int m_BlockSize;
        unsigned int m_BlockSpacing;
        bool m_UseTransformationDam;
        double m_DamDistance;
    };

} // end namespace anima

#include "animaTensorKissingSymmetricBlockMatchingRegistrationMethod.hxx"
