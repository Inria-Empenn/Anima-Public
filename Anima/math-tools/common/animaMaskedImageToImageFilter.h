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
    typedef itk::Image <unsigned char, 3> MaskImageType;
    typedef MaskImageType::RegionType MaskRegionType;
    typedef MaskImageType::Pointer MaskImagePointer;

    /** Set/Get the mask on which to compute estimates. */
    itkSetMacro(ComputationMask, MaskImagePointer)
    itkGetMacro(ComputationMask, MaskImageType *)
    itkSetMacro(ComputationRegion, MaskRegionType)
    itkGetMacro(ComputationRegion, MaskRegionType)

    itkSetMacro(VerboseProgression, bool)

protected:
    MaskedImageToImageFilter()
    {
        m_ComputationMask = NULL;
        m_ComputationRegion.SetSize(0,0);
        m_HighestProcessedSlice = 0;
        m_ProcessedDimension = 0;
        m_VerboseProgression = true;
        m_ProgressReport = 0;
    }

    virtual ~MaskedImageToImageFilter()
    {
        if (m_ProgressReport)
            delete m_ProgressReport;
    }

    virtual void CheckComputationMask();

    void InitializeComputationRegionFromMask();

    virtual void GenerateData();

    static ITK_THREAD_RETURN_TYPE ThreaderMultiSplitCallback(void *arg);
    virtual void ThreadProcessSlices(itk::ThreadIdType threadId);

    virtual void BeforeThreadedGenerateData();

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

    itk::SimpleFastMutexLock m_LockHighestProcessedSlice;
    int m_HighestProcessedSlice;
    unsigned int m_ProcessedDimension;

    itk::ProgressReporter *m_ProgressReport;
    bool m_VerboseProgression;

    MaskImagePointer m_ComputationMask;
    // Optimization of multithread code, compute only on region defined from mask... Uninitialized in constructor.
    MaskRegionType m_ComputationRegion;
};

} //end namespace anima

#include "animaMaskedImageToImageFilter.hxx"
