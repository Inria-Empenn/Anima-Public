#pragma once

#include <iostream>
#include <random>
#include <itkImageToImageFilter.h>
#include <itkImage.h>

#include <time.h>

namespace anima
{

template <class ImageType>
class NoiseGeneratorImageFilter :
public itk::ImageToImageFilter<ImageType,ImageType>
{
public:
    /** Standard class typedefs. */
    typedef NoiseGeneratorImageFilter Self;
    typedef ImageType TInputImage;
    typedef ImageType TOutputImage;
    typedef itk::ImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(NoiseGeneratorImageFilter, ImageToImageFilter)

    typedef typename TOutputImage::PixelType OutputPixelType;
    typedef typename TInputImage::PixelType InputPixelType;

    /** Image typedef support */
    typedef TInputImage  InputImageType;
    typedef TOutputImage OutputImageType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    itkSetMacro(NumberOfReplicates, unsigned int)
    itkGetConstMacro(NumberOfReplicates, unsigned int)
    
    itkSetMacro(StandardDeviation, double)
    itkGetConstMacro(StandardDeviation, double)
    
    itkSetMacro(UseGaussianDistribution, bool)
    itkGetConstMacro(UseGaussianDistribution, bool)

protected:
    NoiseGeneratorImageFilter()
    {
        m_NumberOfReplicates = 1;
        m_StandardDeviation = 1.0;
        m_UseGaussianDistribution = false;
    }

    virtual ~NoiseGeneratorImageFilter()
    {
    }

    void GenerateOutputInformation() ITK_OVERRIDE;
    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(NoiseGeneratorImageFilter);

    unsigned int m_NumberOfReplicates;
    double m_StandardDeviation;
    bool m_UseGaussianDistribution;
    std::vector <std::mt19937> m_Generators;
};

} // end namespace anima

#include "animaNoiseGeneratorImageFilter.hxx"
