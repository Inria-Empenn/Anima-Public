#pragma once

#include <itkImageToImageFilter.h>
#include <animaMCMImage.h>

#include <animaMultiCompartmentModel.h>

namespace anima
{

template <class TPixelType>
class MCMScalarMapsImageFilter :
public itk::ImageToImageFilter< anima::MCMImage <TPixelType, 3>, itk::Image<TPixelType, 3> >
{
public:
    /** Standard class typedefs */
    typedef MCMScalarMapsImageFilter Self;
    typedef anima::MCMImage <TPixelType, 3> InputImageType;
    typedef itk::Image <TPixelType, 3> OutputImageType;
    typedef itk::ImageToImageFilter <InputImageType, OutputImageType > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(MCMScalarMapsImageFilter, itk::ImageToImageFilter)

    typedef typename InputImageType::ConstPointer InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    typedef typename InputImageType::RegionType InputRegionType;
    typedef typename InputImageType::IndexType InputIndexType;
    typedef typename InputImageType::PixelType PixelType;

    // Multi-compartment models typedefs
    typedef anima::MultiCompartmentModel MCModelType;
    typedef MCModelType::Pointer MCModelPointer;

    itkSetMacro(IncludeIsotropicWeights, bool)

protected:
    MCMScalarMapsImageFilter ()
    {
        // Five outputs for now:
        // - free water weight, isotropic restricted weight (sum of IR or Stanisz compartments)
        // - anisotropic weight
        // - FA and MD
        // - others coming soon
        unsigned int numOutputs = 5;
        this->SetNumberOfRequiredOutputs(numOutputs);

        for (unsigned int i = 0;i < numOutputs;++i)
            this->SetNthOutput(i,this->MakeOutput(i));

        m_IncludeIsotropicWeights = false;
    }

    virtual ~MCMScalarMapsImageFilter () {}

    template <class T> bool isZero(const itk::VariableLengthVector <T> &value) const
    {
        for (unsigned int i = 0;i < value.GetNumberOfElements();++i)
        {
            if (value[i] != 0.0)
                return false;
        }

        return true;
    }

    void ThreadedGenerateData(const InputRegionType &region, itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(MCMScalarMapsImageFilter);

    // Use to compute FA and MD measures with or without iso compartments contributions
    bool m_IncludeIsotropicWeights;
};

} // end namespace anima

#include "animaMCMScalarMapsImageFilter.hxx"
