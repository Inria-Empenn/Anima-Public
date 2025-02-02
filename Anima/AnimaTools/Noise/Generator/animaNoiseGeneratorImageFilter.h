#pragma once

#include <iostream>
#include <random>
#include <animaNumberedThreadImageToImageFilter.h>
#include <itkImage.h>

#include <time.h>

namespace anima
{
    template <class ImageType>
    class NoiseGeneratorImageFilter : public anima::NumberedThreadImageToImageFilter<ImageType, itk::Image<double, ImageType::ImageDimension>>
    {
    public:
        /** Standard class typedefs. */
        typedef NoiseGeneratorImageFilter Self;
        typedef ImageType TInputImage;
        typedef itk::Image<double, ImageType::ImageDimension> TOutputImage;
        typedef anima::NumberedThreadImageToImageFilter<TInputImage, TOutputImage> Superclass;
        typedef itk::SmartPointer<Self> Pointer;
        typedef itk::SmartPointer<const Self> ConstPointer;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods) */
        itkTypeMacro(NoiseGeneratorImageFilter, anima::NumberedThreadImageToImageFilter);

        typedef typename TOutputImage::PixelType OutputPixelType;
        typedef typename TInputImage::PixelType InputPixelType;

        /** Image typedef support */
        typedef TInputImage InputImageType;
        typedef TOutputImage OutputImageType;
        typedef typename InputImageType::Pointer InputImagePointer;
        typedef typename OutputImageType::Pointer OutputImagePointer;

        /** Superclass typedefs. */
        typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

        itkSetMacro(NumberOfReplicates, unsigned int);
        itkGetConstMacro(NumberOfReplicates, unsigned int);

        itkSetMacro(NoiseSigma, double);
        itkGetConstMacro(NoiseSigma, double);

        itkSetMacro(UseGaussianDistribution, bool);
        itkGetConstMacro(UseGaussianDistribution, bool);

    protected:
        NoiseGeneratorImageFilter()
        {
            m_NumberOfReplicates = 1;
            m_NoiseSigma = 1.0;
            m_UseGaussianDistribution = false;
        }

        virtual ~NoiseGeneratorImageFilter()
        {
        }

        void GenerateOutputInformation() ITK_OVERRIDE;
        void BeforeThreadedGenerateData() ITK_OVERRIDE;
        void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

    private:
        ITK_DISALLOW_COPY_AND_ASSIGN(NoiseGeneratorImageFilter);

        unsigned int m_NumberOfReplicates;
        double m_NoiseSigma;
        bool m_UseGaussianDistribution;
        std::vector<std::mt19937> m_Generators;
    };

} // end namespace anima

#include "animaNoiseGeneratorImageFilter.hxx"
