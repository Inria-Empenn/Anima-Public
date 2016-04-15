#pragma once

#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>

namespace anima
{

template <class PixelScalarType>
class DTIExtrapolateImageFilter :
    public itk::ImageToImageFilter< itk::VectorImage<PixelScalarType,3>, itk::VectorImage<PixelScalarType,3> >
{
public:
    /** Standard class typedefs. */
    typedef DTIExtrapolateImageFilter<PixelScalarType> Self;
    typedef itk::VectorImage<PixelScalarType,3> TInputImage;
    typedef itk::Image <PixelScalarType,3> OutputB0ImageType;
    typedef itk::VectorImage<PixelScalarType,3> DTIImageType;
    typedef itk::VectorImage<PixelScalarType,3> TOutputImage;
    typedef itk::ImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods) */
    itkTypeMacro(DTIExtrapolateImageFilter, ImageToImageFilter);

    /** Image typedef support */
    typedef TInputImage  InputImageType;
    typedef TOutputImage OutputImageType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    void SetInitialEstimatedB0Image(OutputB0ImageType *b0) {m_InitialEstimatedB0Image = b0;}
    itkGetMacro(EstimatedB0Image, OutputB0ImageType *);

protected:
    DTIExtrapolateImageFilter()
        : Superclass()
    {
        m_InitialEstimatedB0Image = NULL;
    }

    virtual ~DTIExtrapolateImageFilter() {}

    void GenerateOutputInformation() ITK_OVERRIDE;
    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
    DTIExtrapolateImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    typename OutputB0ImageType::Pointer m_InitialEstimatedB0Image;
    typename OutputB0ImageType::Pointer m_EstimatedB0Image;

    static const unsigned int m_NumberOfComponents = 6;
};

} // end namespace anima

#include "animaDTIExtrapolateImageFilter.hxx"
