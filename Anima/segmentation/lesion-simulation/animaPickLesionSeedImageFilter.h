#pragma once

#include <itkImageToImageFilter.h>

namespace anima
{
template <class TInputImage, class TOutputImage>
class PickLesionSeedImageFilter :
public itk::ImageToImageFilter <TInputImage, TOutputImage>
{
public:
    /** Standard class typedefs. */
    typedef PickLesionSeedImageFilter <TInputImage, TOutputImage> Self;
    typedef typename itk::ImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef TInputImage InputImageType;
    typedef TOutputImage OutputImageType;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(PickLesionSeedImageFilter, ImageToImageFilter)

    typedef typename TInputImage::Pointer InputImagePointer;
    typedef typename TInputImage::PixelType InputImagePixel;
    typedef typename TInputImage::IndexType IndexType;
    typedef typename TOutputImage::Pointer OutputImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    itkSetMacro(NumberOfSeeds, unsigned int)
    itkSetMacro(ProximityThreshold, double)

    struct pair_comparator
    {
        bool operator() (const std::pair<IndexType, double> & f, const std::pair<IndexType, double> & s)
        { return (f.second < s.second); }
    };

protected:
    PickLesionSeedImageFilter()
    {
        m_NumberOfSeeds = 1;
        m_ProximityThreshold = 10;
        srand(time(0));
    }

    virtual ~PickLesionSeedImageFilter() {}

    void GenerateData() ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(PickLesionSeedImageFilter);

    unsigned int m_NumberOfSeeds;
    double m_ProximityThreshold;
};

} // end namespace anima

#include "animaPickLesionSeedImageFilter.hxx"
