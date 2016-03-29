#pragma once

#include <iostream>
#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>

#include <random>
#include <time.h>

namespace anima
{

template <unsigned int Dimension>
class NonCentralChiNoiseGeneratorImageFilter :
public itk::ImageToImageFilter< itk::Image <float, Dimension> , itk::Image <float, Dimension> >
{
public:
    /** Standard class typedefs. */
    typedef NonCentralChiNoiseGeneratorImageFilter Self;
    typedef itk::Image <float, Dimension> TInputImage;
    typedef itk::Image <float, Dimension> TOutputImage;
    typedef itk::ImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods) */
    itkTypeMacro(NonCentralChiNoiseGeneratorImageFilter, ImageToImageFilter);

    typedef typename TOutputImage::PixelType OutputPixelType;
    typedef typename TInputImage::PixelType InputPixelType;

    /** Image typedef support */
    typedef TInputImage  InputImageType;
    typedef TOutputImage OutputImageType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    void SetRefVal (double &val)
    {
        m_RefVal = val;
        if (val != 0)
            m_AverageMeanValOnAllVolumes = true;
    }

    itkSetMacro(AverageMeanValOnAllVolumes, bool);
    itkSetMacro(NumberOfCoils, unsigned int);

    itkSetMacro(BackgroundVariance, double);
    itkGetMacro(BackgroundVariance, double);

    itkSetMacro(RelStdDevGaussianNoise, double);
    itkGetMacro(RelStdDevGaussianNoise, double);

protected:
    NonCentralChiNoiseGeneratorImageFilter()
    {
        m_BackgroundVariance = 20;
        m_RelStdDevGaussianNoise = 0.05;
        m_NumberOfCoils = 1;
        m_RefVal = 0;

        m_BackgroundThreshold = 10;
        m_MeanVals.clear();
        m_AverageMeanValOnAllVolumes = false;
    }

    virtual ~NonCentralChiNoiseGeneratorImageFilter()
    {
    }

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId) ITK_OVERRIDE;
    void TreatRegionWithNoiseVariance(const OutputImageRegionType &region, double &variance, itk::ThreadIdType threadId);

private:
    NonCentralChiNoiseGeneratorImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    double m_BackgroundVariance;
    double m_RelStdDevGaussianNoise;

    std::vector <double> m_MeanVals;
    double m_RefVal;

    double m_BackgroundThreshold;
    bool m_AverageMeanValOnAllVolumes;

    //! Number of coils parameter, 1 = Rician noise
    unsigned int m_NumberOfCoils;

    std::vector <std::mt19937> m_Generators;
};
} // end of namespace anima

#include "animaNonCentralChiNoiseGeneratorImageFilter.hxx"
