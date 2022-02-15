#pragma once
#include "animaNumberedThreadImageToImageFilter.h"

namespace anima
{

template <typename TInputImage, typename TOutputImage>
void
NumberedThreadImageToImageFilter <TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
    m_ThreadIdsVector.clear();
    Superclass::BeforeThreadedGenerateData();

    if (m_ComputationRegion.GetSize(0) == 0)
    {
        m_ComputationRegion = this->GetOutput(0)->GetRequestedRegion();
        m_NumberOfPointsToProcess = m_ComputationRegion.GetNumberOfPixels();
    }

    unsigned int maxNumSlices = m_ComputationRegion.GetSize()[0];
    m_ProcessedDimension = 0;
    for (unsigned int i = 1;i < OutputImageRegionType::ImageDimension;++i)
    {
        if (maxNumSlices <= m_ComputationRegion.GetSize()[i])
        {
            maxNumSlices = m_ComputationRegion.GetSize()[i];
            m_ProcessedDimension = i;
        }
    }

    // Since image requested region may now be smaller than the image, fill the outputs with zeros
    typedef typename TOutputImage::PixelType OutputPixelType;

    for (unsigned int i = 0;i < this->GetNumberOfOutputs();++i)
    {
        OutputPixelType zeroPixel;
        this->InitializeZeroPixel(this->GetOutput(i),zeroPixel);
        this->GetOutput(i)->FillBuffer(zeroPixel);
    }
}

template< typename TInputImage, typename TOutputImage >
void
NumberedThreadImageToImageFilter <TInputImage, TOutputImage>
::GenerateData()
{
    this->AllocateOutputs();
    this->BeforeThreadedGenerateData();

    this->ResetMultiThreadingPart();

    ThreadStruct str;
    str.Filter = this;

    this->GetMultiThreader()->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    this->GetMultiThreader()->SetSingleMethod(this->ThreaderMultiSplitCallback, &str);

    this->GetMultiThreader()->SingleMethodExecute();

    this->AfterThreadedGenerateData();
}

template< typename TInputImage, typename TOutputImage >
void
NumberedThreadImageToImageFilter <TInputImage, TOutputImage>
::ResetMultiThreadingPart()
{
    m_HighestProcessedSlice = 0;

    m_NumberOfProcessedPoints = 0;
    this->UpdateProgress(0.0);
}

template< typename TInputImage, typename TOutputImage >
ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION
NumberedThreadImageToImageFilter <TInputImage, TOutputImage>
::ThreaderMultiSplitCallback(void *arg)
{
    ThreadStruct *str;

    str = (ThreadStruct *)( ( (itk::MultiThreaderBase::WorkUnitInfo *)( arg ) )->UserData );

    Self *filterPtr = dynamic_cast <Self *> (str->Filter.GetPointer());
    filterPtr->ThreadProcessSlices();

    return ITK_THREAD_RETURN_DEFAULT_VALUE;
}

template< typename TInputImage, typename TOutputImage >
void
NumberedThreadImageToImageFilter <TInputImage, TOutputImage>
::ThreadProcessSlices()
{
    OutputImageRegionType processedRegion = m_ComputationRegion;
    processedRegion.SetSize(m_ProcessedDimension,1);
    int highestToleratedSliceValue = m_ComputationRegion.GetSize()[m_ProcessedDimension] - 1;

    bool continueLoop = true;
    while (continueLoop)
    {
        m_LockHighestProcessedSlice.lock();

        if (m_HighestProcessedSlice > highestToleratedSliceValue)
        {
            m_LockHighestProcessedSlice.unlock();
            continueLoop = false;
            continue;
        }

        processedRegion.SetIndex(m_ProcessedDimension, m_ComputationRegion.GetIndex()[m_ProcessedDimension] + m_HighestProcessedSlice);
        m_HighestProcessedSlice++;

        m_LockHighestProcessedSlice.unlock();

        this->DynamicThreadedGenerateData(processedRegion);
    }
}

template <typename TInputImage, typename TOutputImage>
void
NumberedThreadImageToImageFilter <TInputImage, TOutputImage>
::IncrementNumberOfProcessedPoints()
{
    m_LockProcessedPoints.lock();

    ++m_NumberOfProcessedPoints;

    double ratio = std::floor(m_NumberOfProcessedPoints * 100.0 / m_NumberOfPointsToProcess) / 100.0;
    ratio = this->progressFixedToFloat(this->progressFloatToFixed(ratio));

    if (ratio != this->GetProgress())
        this->UpdateProgress(ratio);

    m_LockProcessedPoints.unlock();
}

template <typename TInputImage, typename TOutputImage>
unsigned int
NumberedThreadImageToImageFilter <TInputImage, TOutputImage>
::GetSafeThreadId()
{
    m_LockThreadIdNumber.lock();

    unsigned int threadId = 0;
    bool presentInIdsVector = true;
    while (presentInIdsVector)
    {
        presentInIdsVector = false;
        for (unsigned int i = 0;i < m_ThreadIdsVector.size();++i)
        {
            if (m_ThreadIdsVector[i] == threadId)
            {
                presentInIdsVector = true;
                break;
            }
        }

        if (presentInIdsVector)
            threadId++;
    }

    m_ThreadIdsVector.push_back(threadId);

    m_LockThreadIdNumber.unlock();

    return threadId;
}

template <typename TInputImage, typename TOutputImage>
void
NumberedThreadImageToImageFilter <TInputImage, TOutputImage>
::SafeReleaseThreadId(unsigned int threadId)
{
    m_LockThreadIdNumber.lock();

    unsigned int indexThreadId = 0;
    for (unsigned int i = 0;i < m_ThreadIdsVector.size();++i)
    {
        if (threadId == m_ThreadIdsVector[i])
        {
            indexThreadId = i;
            break;
        }
    }

    m_ThreadIdsVector.erase(m_ThreadIdsVector.begin() + indexThreadId);

    m_LockThreadIdNumber.unlock();
}

} // end namespace anima
