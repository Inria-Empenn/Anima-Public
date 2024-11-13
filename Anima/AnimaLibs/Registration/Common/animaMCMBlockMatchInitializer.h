#pragma once
#include <animaBlockMatchInitializer.h>

#include <animaMCMImage.h>
#include <animaMultiCompartmentModel.h>

namespace anima
{

template <class PixelType, unsigned int NDimensions=3>
class MCMBlockMatchingInitializer
        : public anima::BlockMatchingInitializer <PixelType, NDimensions>
{
public:
    typedef anima::BlockMatchingInitializer <PixelType, NDimensions> Superclass;
    typedef MCMBlockMatchingInitializer <PixelType,NDimensions> Self;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    typedef anima::MCMImage <PixelType, NDimensions> MCMImageType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(MCMBlockMatchingInitializer, anima::BlockMatchingInitializer)

    typedef anima::MultiCompartmentModel MCModelType;
    typedef typename MCModelType::Pointer MCModelPointer;

    typedef typename Superclass::BlockGeneratorThreadStruct BlockGeneratorThreadStruct;
    typedef typename Superclass::ImageRegionType ImageRegionType;

    struct MCMBlockGeneratorThreadStruct : public BlockGeneratorThreadStruct
    {
        std::vector < std::vector <MCModelPointer> > ref_models;
    };

    void AddReferenceImage(itk::ImageBase <NDimensions> *refImage) ITK_OVERRIDE;

protected:
    MCMBlockMatchingInitializer() : Superclass()
    {
    }

    virtual ~MCMBlockMatchingInitializer() {}

    void InitializeThreading(unsigned int maskIndex, BlockGeneratorThreadStruct *&workStr) ITK_OVERRIDE;
    bool CheckOrientedModelVariance(unsigned int imageIndex, ImageRegionType &region, double &blockVariance,
                                    BlockGeneratorThreadStruct *workStr, unsigned int threadId) ITK_OVERRIDE;

private:
    MCMBlockMatchingInitializer(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    std::vector <MCModelPointer> m_ReferenceModels;
};

} // end of namespace anima

#include "animaMCMBlockMatchInitializer.hxx"
