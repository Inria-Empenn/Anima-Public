#pragma once

#include <itkImageToImageFilter.h>
#include <animaMCMImage.h>

#include <animaMultiCompartmentModel.h>
#include <animaMCMWeightedAverager.h>

namespace anima
{

template <class TPixelType>
class MCMAverageImagesImageFilter :
public itk::ImageToImageFilter< anima::MCMImage <TPixelType, 3>, anima::MCMImage<TPixelType, 3> >
{
public:
    /** Standard class type def */

    typedef MCMAverageImagesImageFilter Self;
    typedef anima::MCMImage <TPixelType, 3> InputImageType;
    typedef anima::MCMImage <TPixelType, 3> OutputImageType;
    typedef itk::ImageToImageFilter <InputImageType, OutputImageType > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(MCMAverageImagesImageFilter, ImageToImageFilter)

    typedef typename InputImageType::ConstPointer InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    typedef typename InputImageType::RegionType InputRegionType;
    typedef typename InputImageType::IndexType InputIndexType;
    typedef typename InputImageType::PixelType PixelType;

    // Multi-compartment models typedefs
    typedef anima::MultiCompartmentModel MCModelType;
    typedef MCModelType::Pointer MCModelPointer;

    typedef anima::MCMWeightedAverager MCMAveragerType;
    typedef typename MCMAveragerType::Pointer MCMAveragerPointer;

    void SetReferenceOutputModel(MCModelPointer &model);
    MCModelType *GetReferenceOutputModel() {return m_ReferenceOutputModel;}

protected:
    MCMAverageImagesImageFilter ()
    {
    }

    virtual ~MCMAverageImagesImageFilter () {}

    bool isZero(const itk::VariableLengthVector <double> &value) const
    {
        for (unsigned int i = 0;i < value.GetNumberOfElements();++i)
        {
            if (value[i] != 0)
                return false;
        }

        return true;
    }

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void DynamicThreadedGenerateData(const InputRegionType &region) ITK_OVERRIDE;
    virtual MCMAveragerPointer CreateAverager();

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(MCMAverageImagesImageFilter);

    std::vector <MCModelPointer> m_ReferenceInputModels;
    MCModelPointer m_ReferenceOutputModel;
};

} // end namespace anima

#include "animaMCMAverageImagesImageFilter.hxx"
