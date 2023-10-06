#pragma once

#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>

#include <animaMCMImage.h>
#include <animaMultiCompartmentModel.h>

namespace anima
{

class DDITestAveragingOnRealValueImageFilter :
public itk::ImageToImageFilter< anima::MCMImage <double, 3>, anima::MCMImage<double, 3> >
{
public:

    /** Standard class type def */
    typedef DDITestAveragingOnRealValueImageFilter Self;
    typedef anima::MCMImage <double, 3> InputImageType;
    typedef anima::MCMImage <double, 3> OutputImageType;
    typedef itk::ImageToImageFilter <InputImageType, OutputImageType > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(DDITestAveragingOnRealValueImageFilter, ImageToImageFilter)

    typedef InputImageType::ConstPointer InputImagePointer;
    typedef OutputImageType::Pointer     OutputImagePointer;

    typedef InputImageType::RegionType InputRegionType;
    typedef InputImageType::IndexType  InputIndexType;
    typedef InputImageType::PixelType  PixelType;

    // Multi-compartment models typedefs
    typedef anima::MultiCompartmentModel MCModelType;
    typedef MCModelType::Pointer MCModelPointer;

    itkSetMacro(Step, int)
    itkSetMacro(Method, int)

    void SetReferenceOutputModel(MCModelPointer &model);

protected:
    DDITestAveragingOnRealValueImageFilter ()
    {
        m_Step = 2;
    }
    virtual ~DDITestAveragingOnRealValueImageFilter () {}

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void DynamicThreadedGenerateData(const InputRegionType &region) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(DDITestAveragingOnRealValueImageFilter);

    int m_Step;
    int m_Method;

    MCModelPointer m_ReferenceInputModel;
    MCModelPointer m_ReferenceOutputModel;
};

}
