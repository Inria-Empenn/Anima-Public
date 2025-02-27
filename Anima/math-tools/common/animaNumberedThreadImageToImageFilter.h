#pragma once

#include <itkImageToImageFilter.h>
#include <mutex>
#include <itkVariableLengthVector.h>

namespace anima
{

/**
 * @brief Implements a class to handle thread number in a dynamic way for multithreaded methods needing
 * thread numbering even for dynamic threading
 */
template <typename TInputImage, typename TOutputImage>
class NumberedThreadImageToImageFilter :
        public itk::ImageToImageFilter < TInputImage, TOutputImage >
{
public:
    /** Standard class typedefs. */
    typedef NumberedThreadImageToImageFilter Self;
    typedef itk::ImageToImageFilter <TInputImage, TOutputImage> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;
    typedef typename itk::ImageSource<TOutputImage>::ThreadStruct ThreadStruct;

    static constexpr unsigned int OutputImageDimension = TOutputImage::ImageDimension;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(NumberedThreadImageToImageFilter, itk::ImageToImageFilter)

    itkSetMacro(NumberOfPointsToProcess, unsigned int)

    itkSetMacro(ComputationRegion, OutputImageRegionType)
    itkGetMacro(ComputationRegion, OutputImageRegionType)

protected:
    NumberedThreadImageToImageFilter()
    {
        m_NumberOfProcessedPoints = 0;
        m_ComputationRegion.SetSize(0,0);
        m_HighestProcessedSlice = 0;
        m_ProcessedDimension = 0;
    }

    virtual ~NumberedThreadImageToImageFilter() {}

    virtual void GenerateData() ITK_OVERRIDE;
    virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;

    static ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION ThreaderMultiSplitCallback(void *arg);
    virtual void ThreadProcessSlices();

    unsigned int GetSafeThreadId();
    void SafeReleaseThreadId(unsigned int threadId);

    void IncrementNumberOfProcessedPoints();
    void ResetMultiThreadingPart();

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
    ITK_DISALLOW_COPY_AND_ASSIGN(NumberedThreadImageToImageFilter);

    std::mutex m_LockThreadIdNumber;
    std::vector <unsigned int> m_ThreadIdsVector;

    std::mutex m_LockProcessedPoints;
    unsigned int m_NumberOfProcessedPoints;
    unsigned int m_NumberOfPointsToProcess;

    std::mutex m_LockHighestProcessedSlice;
    int m_HighestProcessedSlice;
    unsigned int m_ProcessedDimension;

    // Optimization of multithread code, compute only on region defined from mask... Uninitialized in constructor.
    OutputImageRegionType m_ComputationRegion;
};

} //end namespace anima

#include "animaNumberedThreadImageToImageFilter.hxx"
