#pragma once

#include <iostream>
#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>
#include <vector>

namespace anima
{

template <typename TInputPixelType>
class GeneralizedFAImageFilter :
public itk::ImageToImageFilter< itk::VectorImage<TInputPixelType, 3> , itk::Image <TInputPixelType, 3> >
{
public:
    /** Standard class typedefs. */
    typedef GeneralizedFAImageFilter Self;
    typedef itk::VectorImage <TInputPixelType, 3> TInputImage;
    typedef itk::Image <TInputPixelType, 3> TOutputImage;
    typedef itk::ImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods) */
    itkTypeMacro(GeneralizedFAImageFilter, ImageToImageFilter);

    typedef typename TInputImage::Pointer InputImagePointer;
    typedef typename TInputImage::PixelType InputImagePixel;
    typedef typename TOutputImage::Pointer OutputImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

protected:
    GeneralizedFAImageFilter()
    {
    }

    virtual ~GeneralizedFAImageFilter() {}

    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId) ITK_OVERRIDE;

    bool isZero(InputImagePixel &testVal)
    {
        bool resVal = true;
        for (unsigned int i = 0;i < testVal.GetSize();++i)
        {
            if (testVal[i] != 0)
            {
                resVal = false;
                break;
            }
        }

        return resVal;
    }

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(GeneralizedFAImageFilter);
};
	
} // end of namespace anima

#include "animaGeneralizedFAImageFilter.hxx"
