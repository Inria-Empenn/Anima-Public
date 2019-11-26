#pragma once

#include <iostream>
#include <animaNumberedThreadImageToImageFilter.h>
#include <itkNumericTraits.h>

namespace anima
{

template <typename TInputImage, typename TOutputImage>
class MaskedImageToImageFilter :
        public anima::NumberedThreadImageToImageFilter <TInputImage, TOutputImage>
{
public:
    /** Standard class typedefs. */
    typedef MaskedImageToImageFilter Self;
    typedef anima::NumberedThreadImageToImageFilter <TInputImage, TOutputImage> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(MaskedImageToImageFilter, anima::NumberedThreadImageToImageFilter)

    /** Superclass typedefs. */
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
    typedef typename itk::ImageSource<TOutputImage>::ThreadStruct ThreadStruct;

    /** Mask typedefs */
    typedef itk::Image <unsigned char,TInputImage::ImageDimension> MaskImageType;
    typedef typename MaskImageType::RegionType MaskRegionType;
    typedef typename MaskImageType::Pointer MaskImagePointer;
    typedef typename MaskImageType::IndexType MaskIndexType;
    typedef typename MaskImageType::SizeType MaskSizeType;

    /** Set/Get the mask on which to compute estimates. */
    itkSetMacro(ComputationMask, MaskImagePointer)
    itkGetMacro(ComputationMask, MaskImageType *)

protected:
    MaskedImageToImageFilter()
    {
        m_ComputationMask = ITK_NULLPTR;
    }

    virtual ~MaskedImageToImageFilter() {}

    virtual void CheckComputationMask();
    void InitializeComputationRegionFromMask();
    virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(MaskedImageToImageFilter);

    MaskImagePointer m_ComputationMask;
};

} //end namespace anima

#include "animaMaskedImageToImageFilter.hxx"
