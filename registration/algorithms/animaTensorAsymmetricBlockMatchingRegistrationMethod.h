#pragma once

#include <animaTensorBlockMatchingBaseRegistrationMethod.h>

namespace anima
{

    template <typename TInputImage>
    class TensorAsymmetricBlockMatchingRegistrationMethod : public TensorBlockMatchingBaseRegistrationMethod <TInputImage>
    {
    public:
        /** Standard class typedefs. */
        typedef TensorAsymmetricBlockMatchingRegistrationMethod Self;
        typedef TensorBlockMatchingBaseRegistrationMethod <TInputImage> Superclass;
        typedef itk::SmartPointer<Self> Pointer;
        typedef itk::SmartPointer<const Self> ConstPointer;
        
        typedef typename Superclass::InputImageType InputImageType;
        typedef typename Superclass::TransformType TransformType;
        typedef typename Superclass::TransformPointer TransformPointer;
        typedef typename Superclass::ThreadedMatchData ThreadedMatchData;
        
        /** Method for creation through the object factory. */
        itkNewMacro(Self);
        
        /** Run-time type information (and related methods). */
        itkTypeMacro(TensorAsymmetricBlockMatchingRegistrationMethod, BlockMatchingBaseRegistrationMethod);
        
    protected:
        TensorAsymmetricBlockMatchingRegistrationMethod();
        virtual ~TensorAsymmetricBlockMatchingRegistrationMethod() {}
        
        virtual void PerformOneIteration(InputImageType *refImage, InputImageType *movingImage,
                                         TransformPointer &addOn);
        
    private:
        TensorAsymmetricBlockMatchingRegistrationMethod(const Self&); //purposely not implemented
        void operator=(const Self&); //purposely not implemented
    };
    
} // end namespace anima

#include "animaTensorAsymmetricBlockMatchingRegistrationMethod.hxx"
