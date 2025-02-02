#pragma once

#include <itkImageToImageFilter.h>
#include <itkImage.h>

namespace anima
{

template <typename TInputImage, typename TOutputImage>
class PseudoResidualsNoiseImageFilter :
        public itk::ImageToImageFilter< TInputImage, TOutputImage >
{
public:
    /** Standard class typedefs. */
    typedef PseudoResidualsNoiseImageFilter Self;
    typedef itk::ImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(PseudoResidualsNoiseImageFilter, ImageToImageFilter)

    /** Image typedef support */
    typedef TInputImage  InputImageType;
    typedef TOutputImage OutputImageType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename InputImageType::IndexType InputImageIndexType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    itkSetMacro(PatchHalfSize, unsigned int)

protected:
    PseudoResidualsNoiseImageFilter()
    : Superclass()
    {
        m_PatchHalfSize = 1;
    }

    virtual ~PseudoResidualsNoiseImageFilter() {}

    void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(PseudoResidualsNoiseImageFilter);

    unsigned int m_PatchHalfSize;
};

} // end namespace anima

#include "animaPseudoResidualsNoiseImageFilter.hxx"
