#pragma once

#include <itkImageToImageFilter.h>

namespace anima
{
template <typename TInputPixelType>
class Bootstrap4DVolumeImageFilter :
        public itk::ImageToImageFilter< itk::Image<TInputPixelType, 4> , itk::Image <TInputPixelType, 4> >
{
public:
    /** Standard class typedefs. */
    typedef Bootstrap4DVolumeImageFilter Self;
    typedef itk::Image <TInputPixelType, 4> TInputImage;
    typedef itk::Image <TInputPixelType, 4> TOutputImage;
    typedef itk::ImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(Bootstrap4DVolumeImageFilter, ImageToImageFilter)

    typedef typename TInputImage::Pointer InputImagePointer;
    typedef typename TInputImage::PixelType InputImagePixel;
    typedef typename TOutputImage::Pointer OutputImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    void SetNthInformationFile(unsigned int i, std::string fileName);

    std::vector <std::string> &GetOutputInformation() {return m_OutputInformation;}

protected:
    Bootstrap4DVolumeImageFilter()
    {
        m_ImageInformations.clear();
        m_OutputInformation.clear();
    }

    virtual ~Bootstrap4DVolumeImageFilter() {}

    void GenerateData() ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(Bootstrap4DVolumeImageFilter);

    std::vector < std::vector <std::string> > m_ImageInformations;
    std::vector <std::string> m_OutputInformation;
};

} // end of namespace anima

#include "animaBootstrap4DVolumeImageFilter.hxx"
