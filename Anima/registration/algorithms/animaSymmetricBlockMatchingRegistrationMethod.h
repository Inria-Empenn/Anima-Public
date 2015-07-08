#pragma once

#include <animaBlockMatchingBaseRegistrationMethod.h>

namespace anima
{

    template <typename TInputImage>
    class SymmetricBlockMatchingRegistrationMethod : public BlockMatchingBaseRegistrationMethod <TInputImage>
    {
    public:
        /** Standard class typedefs. */
        typedef SymmetricBlockMatchingRegistrationMethod Self;
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
        
        /** Method for creation through the object factory. */
        itkNewMacro(Self);
        
        /** Run-time type information (and related methods). */
        itkTypeMacro(SymmetricBlockMatchingRegistrationMethod, BlockMatchingBaseRegistrationMethod);
        
        /** Set / Get Block Origins */
        void SetFixedBlockRegions(std::vector <ImageRegionType>& origins) {m_FixedBlockRegions = origins;}
        std::vector <ImageRegionType>& GetFixedBlockRegions() {return m_FixedBlockRegions;}
        void SetMovingBlockRegions(std::vector <ImageRegionType>& origins) {m_MovingBlockRegions = origins;}
        std::vector <ImageRegionType>& GetMovingBlockRegions() {return m_MovingBlockRegions;}
        
        /** Set/Get the dam indexes. */
        void SetFixedDamIndexes(std::vector <ImageIndexType>& damInd) {m_FixedDamIndexes = damInd;}
        std::vector <ImageIndexType>& GetFixedDamIndexes() {return m_FixedDamIndexes;}
        
        void SetMovingDamIndexes(std::vector <ImageIndexType>& damInd) {m_MovingDamIndexes = damInd;}
        std::vector <ImageIndexType>& GetMovingDamIndexes() {return m_MovingDamIndexes;}

        
    protected:
        SymmetricBlockMatchingRegistrationMethod();
        virtual ~SymmetricBlockMatchingRegistrationMethod() {}
        
        virtual void PerformOneIteration(InputImageType *refImage, InputImageType *movingImage, TransformPointer &addOn);
        
    private:
        SymmetricBlockMatchingRegistrationMethod(const Self&); //purposely not implemented
        void operator=(const Self&); //purposely not implemented
        
        std::vector<ImageRegionType> m_FixedBlockRegions, m_MovingBlockRegions;
        std::vector <ImageIndexType> m_FixedDamIndexes, m_MovingDamIndexes;
    };
    
} // end namespace anima

#include "animaSymmetricBlockMatchingRegistrationMethod.hxx"
