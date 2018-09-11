#pragma once

#include <itkImageToImageFilter.h>
#include <itkFastMutexLock.h>

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

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(NumberedThreadImageToImageFilter, itk::ImageToImageFilter)

protected:
    NumberedThreadImageToImageFilter() {}

    virtual ~NumberedThreadImageToImageFilter() {}

    virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;

    unsigned int GetSafeThreadId();
    void SafeReleaseThreadId(unsigned int threadId);

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(NumberedThreadImageToImageFilter);

    itk::SimpleFastMutexLock m_LockThreadIdNumber;
    std::vector <unsigned int> m_ThreadIdsVector;
};

} //end namespace anima

#include "animaNumberedThreadImageToImageFilter.hxx"
