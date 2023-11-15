#pragma once

#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>

namespace anima
{

    template <typename TInputPixelType>
    class SampleImageFromDistributionImageFilter : public itk::ImageToImageFilter<itk::VectorImage<TInputPixelType, 3>, itk::VectorImage<TInputPixelType, 3>>
    {
    public:
        /** Standard class typedefs. */
        typedef SampleImageFromDistributionImageFilter Self;
        typedef itk::VectorImage<TInputPixelType, 3> TInputImage;
        typedef itk::VectorImage<TInputPixelType, 3> TOutputImage;
        typedef itk::ImageToImageFilter<TInputImage, TOutputImage> Superclass;
        typedef itk::SmartPointer<Self> Pointer;
        typedef itk::SmartPointer<const Self> ConstPointer;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods) */
        itkTypeMacro(SampleImageFromDistributionImageFilter, ImageToImageFilter);

        typedef typename TInputImage::Pointer InputImagePointer;
        typedef typename TInputImage::PixelType InputImagePixel;
        typedef typename TOutputImage::Pointer OutputImagePointer;

        /** Superclass typedefs. */
        typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    protected:
        SampleImageFromDistributionImageFilter()
        {
            m_VectorSize = 3;
        }

        virtual ~SampleImageFromDistributionImageFilter() {}

        void GenerateOutputInformation() ITK_OVERRIDE;
        void BeforeThreadedGenerateData() ITK_OVERRIDE;
        void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

    private:
        ITK_DISALLOW_COPY_AND_ASSIGN(SampleImageFromDistributionImageFilter);

        InputImagePixel m_BaseDistributionSample;
        unsigned int m_VectorSize;
    };

} // end namespace anima

#include "animaSampleImageFromDistributionImageFilter.hxx"
