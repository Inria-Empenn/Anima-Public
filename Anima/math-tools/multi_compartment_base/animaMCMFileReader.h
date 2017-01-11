#pragma once

#include <animaMultiCompartmentModel.h>
#include <animaMCMImage.h>
#include <itkImageIOBase.h>

#include <string>
#include <animaBaseCompartment.h>

#include <AnimaMCMBaseExport.h>

namespace anima
{

// For distinguishing image internal pixel type, used for medInria plugins
itk::ImageIOBase::IOComponentType ANIMAMCMBASE_EXPORT GetMCMComponentType(std::string fileName);

template <class PixelType, unsigned int ImageDimension>
class MCMFileReader
{
public:
    typedef anima::MultiCompartmentModel ModelType;
    typedef ModelType::Pointer ModelPointer;
    typedef anima::MCMImage <PixelType, ImageDimension> OutputImageType;
    typedef typename OutputImageType::Pointer OutputImagePointer;
    typedef itk::VectorImage <PixelType, ImageDimension> BaseInputImageType;
    typedef typename BaseInputImageType::Pointer BaseInputImagePointer;

    MCMFileReader();
    ~MCMFileReader();

    OutputImagePointer &GetModelVectorImage() {return m_OutputImage;}
    void SetFileName(std::string fileName) {m_FileName = fileName;}

    void Update();
    virtual anima::BaseCompartment::Pointer CreateCompartmentForType(std::string &compartmentType);

private:
    OutputImagePointer m_OutputImage;
    std::string m_FileName;
};

} // end namespace anima

#include "animaMCMFileReader.hxx"
