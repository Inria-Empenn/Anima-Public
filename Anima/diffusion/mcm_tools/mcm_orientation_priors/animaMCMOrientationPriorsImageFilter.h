#pragma once

#include <animaMCMImage.h>
#include <animaMultiCompartmentModel.h>

#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>

namespace anima
{

template <class TPixelType>
class MCMOrientationPriorsImageFilter :
public itk::ImageToImageFilter< anima::MCMImage <TPixelType, 3>, itk::VectorImage<TPixelType, 3> >
{
public:
    /** Standard class type def */

    using Self = MCMOrientationPriorsImageFilter;
    using InputImageType = anima::MCMImage <TPixelType, 3>;
    using OutputImageType = itk::VectorImage <TPixelType, 3>;
    using Superclass = itk::ImageToImageFilter <InputImageType, OutputImageType >;
    using Pointer = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer<const Self>;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(MCMOrientationPriorsImageFilter, ImageToImageFilter)

    using InputImagePointer = typename InputImageType::ConstPointer;
    using InputRegionType = typename InputImageType::RegionType;

    using OutputImagePointer = typename OutputImageType::Pointer;
    using OutputPixelType = typename OutputImageType::PixelType;

    using MaskImageType = itk::Image <unsigned char, 3>;
    using MaskImagePointer = MaskImageType::Pointer;

    // Multi-compartment models typedefs
    using MCModelType = anima::MultiCompartmentModel;
    using MCModelPointer = MCModelType::Pointer;

    void AddMaskImage(MaskImageType *maskImage) {m_MaskImages.push_back(maskImage);}
    unsigned int GetNumberOfAnisotropicCompartments() {return m_NumberOfAnisotropicCompartments;}

protected:
    MCMOrientationPriorsImageFilter()
    {
        m_ReferenceInputModels.clear();
        m_MaskImages.clear();
        m_NumberOfAnisotropicCompartments = 0;
    }

    virtual ~MCMOrientationPriorsImageFilter() {}

    bool isZero(const itk::VariableLengthVector <double> &value) const
    {
        for (unsigned int i = 0;i < value.GetNumberOfElements();++i)
        {
            if (value[i] != 0)
                return false;
        }

        return true;
    }

    void GenerateOutputInformation() ITK_OVERRIDE;
    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void DynamicThreadedGenerateData(const InputRegionType &region) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(MCMOrientationPriorsImageFilter);

    std::vector <MCModelPointer> m_ReferenceInputModels;
    std::vector <MaskImagePointer> m_MaskImages;
    unsigned int m_NumberOfAnisotropicCompartments;
};

} // end namespace anima

#include "animaMCMOrientationPriorsImageFilter.hxx"
