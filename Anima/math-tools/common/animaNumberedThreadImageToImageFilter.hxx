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

    m_NumberOfProcessedPoints = 0;
    m_NumberOfPointsToProcess = this->GetOutput(0)->GetRequestedRegion().GetNumberOfPixels();
    this->UpdateProgress(0.0);
}

template <typename TInputImage, typename TOutputImage>
void
NumberedThreadImageToImageFilter <TInputImage, TOutputImage>
::IncrementNumberOfProcessedPoints()
{
    m_LockProcessedPoints.Lock();

    ++m_NumberOfProcessedPoints;

    double ratio = m_NumberOfProcessedPoints / m_NumberOfPointsToProcess;
    this->UpdateProgress(ratio);

    m_LockProcessedPoints.Unlock();
}

template <typename TInputImage, typename TOutputImage>
unsigned int
NumberedThreadImageToImageFilter <TInputImage, TOutputImage>
::GetSafeThreadId()
{
    m_LockThreadIdNumber.Lock();

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

    m_LockThreadIdNumber.Unlock();

    return threadId;
}

template <typename TInputImage, typename TOutputImage>
void
NumberedThreadImageToImageFilter <TInputImage, TOutputImage>
::SafeReleaseThreadId(unsigned int threadId)
{
    m_LockThreadIdNumber.Lock();

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

    m_LockThreadIdNumber.Unlock();
}

} // end namespace anima
