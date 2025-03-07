#pragma once

#include <animaMultiCompartmentModel.h>
#include <animaMCMImage.h>

namespace anima
{

template <class PixelType, unsigned int ImageDimension>
class MCMFileWriter
{
public:
    typedef anima::MultiCompartmentModel ModelType;
    typedef ModelType::Pointer ModelPointer;
    typedef anima::MCMImage <PixelType, ImageDimension> InputImageType;
    typedef typename InputImageType::Pointer InputImagePointer;

    typedef itk::VectorImage <PixelType, ImageDimension> BaseOutputImageType;
    typedef typename BaseOutputImageType::Pointer BaseOutputImagePointer;

    MCMFileWriter();
    ~MCMFileWriter();

    void SetInputImage(InputImageType *input) {m_InputImage = input;}
    void SetFileName(std::string fileName);

    void Update();

private:
    InputImagePointer m_InputImage;
    std::string m_FileName;
};

} // end namespace anima

#include "animaMCMFileWriter.hxx"
