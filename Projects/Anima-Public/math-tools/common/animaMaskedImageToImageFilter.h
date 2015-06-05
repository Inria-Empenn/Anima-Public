#pragma once

#include <iostream>
#include <itkImageToImageFilter.h>
#include <animaImageRegionSplitterMask.h>
#include <itkNumericTraits.h>
#include <itkVariableLengthVector.h>

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
    itkNewMacro(Self);

    /** Run-time type information (and related methods) */
    itkTypeMacro(MaskedImageToImageFilter, itk::ImageToImageFilter);

    /** Superclass typedefs. */
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    /** Mask typedefs */
    typedef itk::Image <unsigned char, 3> MaskImageType;
    typedef MaskImageType::RegionType MaskRegionType;
    typedef MaskImageType::Pointer MaskImagePointer;

    typedef anima::ImageRegionSplitterMask <MaskImageType> RegionSplitterType;
    typedef typename RegionSplitterType::Pointer RegionSplitterPointer;

    /** Set/Get the mask on which to compute STAPLE estimate. */
    itkSetMacro(ComputationMask, MaskImagePointer);
    itkGetMacro(ComputationMask, MaskImageType *);
    itkSetMacro(ComputationRegion, MaskRegionType);
    itkGetMacro(ComputationRegion, MaskRegionType);

    unsigned int GetActualNumberOfThreads()
    {
        return m_RegionSplitter->GetNumberOfThreads();
    }

protected:
    MaskedImageToImageFilter()
    {
        m_ComputationMask = NULL;

        m_ComputationRegion.SetSize(0,0);
        m_RegionSplitter = RegionSplitterType::New();
    }

    virtual ~MaskedImageToImageFilter() {}

    virtual void CheckComputationMask();

    void InitializeComputationRegionFromMask();

    virtual const itk::ImageRegionSplitterBase* GetImageRegionSplitter(void) const
    {
        return m_RegionSplitter;
    }

    virtual void BeforeThreadedGenerateData(void);

    virtual void PrintSelf(std::ostream& os, itk::Indent indent) const
    {
        Superclass::PrintSelf(os,indent);
    }

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
    MaskedImageToImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    //For overriding the split image region with something more intelligent taking into account the mask
    RegionSplitterPointer m_RegionSplitter;

    MaskImagePointer m_ComputationMask;
    // Optimization of multithread code, compute only on region defined from mask... Uninitialized in constructor.
    MaskRegionType m_ComputationRegion;
};

} //end namespace anima

#include "animaMaskedImageToImageFilter.hxx"
