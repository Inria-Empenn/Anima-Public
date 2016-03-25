#pragma once

#include <itkImageToImageFilter.h>

namespace anima
{

template <class PixelScalarType>
class FDRCorrectImageFilter :
public itk::ImageToImageFilter< itk::Image<PixelScalarType, 3> , itk::Image <unsigned char, 3> >
{
public:
    /** Standard class typedefs. */
    typedef FDRCorrectImageFilter<PixelScalarType> Self;
    typedef itk::Image <PixelScalarType, 3> TInputImage;
    typedef itk::Image <unsigned char, 3> TOutputImage;
    typedef itk::ImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(FDRCorrectImageFilter, ImageToImageFilter)

    typedef typename TInputImage::Pointer InputImagePointer;
    typedef typename TInputImage::PixelType InputImagePixel;
    typedef typename TOutputImage::Pointer OutputImagePointer;

    typedef itk::Image <unsigned char, 3> MaskImageType;
    typedef typename MaskImageType::Pointer MaskImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    itkSetMacro(MaskImage, MaskImagePointer)
    itkSetMacro(QValue, double)
    itkSetMacro(BYCorrection, bool)

protected:
    FDRCorrectImageFilter()
    {
        m_MaskImage = NULL;
        m_QValue = 0.05;
        m_BYCorrection = false;
    }

    virtual ~FDRCorrectImageFilter() {}

    void GenerateData();

    void PrintSelf(std::ostream& os, itk::Indent indent) const
    {
        Superclass::PrintSelf(os,indent);
    }

private:
    FDRCorrectImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    void CreateFullMask();

    MaskImagePointer m_MaskImage;
    double m_QValue;
    bool m_BYCorrection;
};

} // end of namespace anima

#include "animaFDRCorrectImageFilter.hxx"
