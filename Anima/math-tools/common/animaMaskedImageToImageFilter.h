#pragma once

#include <iostream>
#include <itkImageToImageFilter.h>
#include <itkNumericTraits.h>
#include <itkVariableLengthVector.h>
#include <itkFastMutexLock.h>
#include <itkProgressReporter.h>

namespace anima
{

template <typename TInputImage, typename TOutputImage>
class MaskedImageToImageFilter :
        public itk::ImageToImageFilter < TInputImage, TOutputImage >
{
public:
    /** Standard class typedefs. */
    typedef MaskedImageToImageFilter Self;
    typedef itk::ImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(MaskedImageToImageFilter, itk::ImageToImageFilter)

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
    virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;

    //! Utility function to initialize output images pixel to zero for vector images
    template <typename ScalarRealType>
    void
    InitializeZeroPixel(TOutputImage *image, itk::VariableLengthVector <ScalarRealType> &zeroPixel)
    {
        unsigned int vectorSize = image->GetVectorLength();
        zeroPixel = itk::VariableLengthVector <ScalarRealType> (vectorSize);
        zeroPixel.Fill(0.0);
    }

    //! Utility function to initialize output images pixel to zero for all images except vector images
    template <typename PixelType>
    void
    InitializeZeroPixel(TOutputImage *itkNotUsed(image), PixelType &zeroPixel)
    {
        zeroPixel = itk::NumericTraits <PixelType>::ZeroValue();
    }

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(MaskedImageToImageFilter);

    MaskImagePointer m_ComputationMask;
};

} //end namespace anima

#include "animaMaskedImageToImageFilter.hxx"
