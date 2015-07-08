#pragma once

#include <animaBlockMatchingBaseRegistrationMethod.h>

namespace anima
{

    template <typename TInputImage>
    class AsymmetricBlockMatchingRegistrationMethod : public BlockMatchingBaseRegistrationMethod <TInputImage>
    {
    public:
        /** Standard class typedefs. */
        typedef AsymmetricBlockMatchingRegistrationMethod Self;
        typedef BlockMatchingBaseRegistrationMethod <TInputImage> Superclass;
        typedef itk::SmartPointer<Self> Pointer;
        typedef itk::SmartPointer<const Self> ConstPointer;
        
        typedef typename Superclass::InputImageType InputImageType;
        typedef typename Superclass::TransformType TransformType;
        typedef typename Superclass::TransformPointer TransformPointer;
        typedef typename Superclass::ThreadedMatchData ThreadedMatchData;
        
        /** Method for creation through the object factory. */
        itkNewMacro(Self);
        
        /** Run-time type information (and related methods). */
        itkTypeMacro(AsymmetricBlockMatchingRegistrationMethod, BlockMatchingBaseRegistrationMethod);
        
    protected:
        AsymmetricBlockMatchingRegistrationMethod();
        virtual ~AsymmetricBlockMatchingRegistrationMethod() {}
        
        virtual void PerformOneIteration(InputImageType *refImage, InputImageType *movingImage,
                                         TransformPointer &addOn);
        
    private:
        AsymmetricBlockMatchingRegistrationMethod(const Self&); //purposely not implemented
        void operator=(const Self&); //purposely not implemented
    };
    
} // end namespace anima

#include "animaAsymmetricBlockMatchingRegistrationMethod.hxx"
