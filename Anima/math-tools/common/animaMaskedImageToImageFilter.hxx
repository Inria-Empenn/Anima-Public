#pragma once
#include "animaMaskedImageToImageFilter.h"
#include <itkMultiThreader.h>

#include <itkImageRegionConstIterator.h>

namespace anima
{
template< typename TInputImage, typename TOutputImage >
void
MaskedImageToImageFilter < TInputImage, TOutputImage >
::CheckComputationMask()
{
    if (!m_ComputationMask)
    {
        m_ComputationMask = MaskImageType::New();
        m_ComputationMask->Initialize();

        m_ComputationMask->SetRegions (this->GetOutput(0)->GetLargestPossibleRegion());
        m_ComputationMask->SetSpacing (this->GetOutput(0)->GetSpacing());
        m_ComputationMask->SetOrigin (this->GetOutput(0)->GetOrigin());
        m_ComputationMask->SetDirection (this->GetOutput(0)->GetDirection());

        m_ComputationMask->Allocate();

        m_ComputationMask->FillBuffer(1);
    }
}

template< typename TInputImage, typename TOutputImage >
void
MaskedImageToImageFilter < TInputImage, TOutputImage >
::InitializeComputationRegionFromMask()
{
    this->CheckComputationMask();

    typedef itk::ImageRegionConstIterator< MaskImageType > MaskRegionIteratorType;

    MaskRegionIteratorType maskItr(m_ComputationMask,m_ComputationMask->GetLargestPossibleRegion());
    maskItr.GoToBegin();

    MaskIndexType minPos, maxPos;

    for (unsigned int i = 0;i < m_ComputationMask->GetImageDimension();++i)
    {
        minPos[i] = m_ComputationMask->GetLargestPossibleRegion().GetIndex()[i] + m_ComputationMask->GetLargestPossibleRegion().GetSize()[i];
        maxPos[i] = 0;
    }

    while (!maskItr.IsAtEnd())
    {
        if (maskItr.Get() != 0)
        {
            MaskIndexType tmpInd = maskItr.GetIndex();

            for (unsigned int i = 0;i < m_ComputationMask->GetImageDimension();++i)
            {
                if (minPos[i] > tmpInd[i])
                    minPos[i] = tmpInd[i];

                if (maxPos[i] < tmpInd[i])
                    maxPos[i] = tmpInd[i];
            }
        }

        ++maskItr;
    }

    m_ComputationRegion.SetIndex(minPos);

    MaskSizeType tmpSize;
    for (unsigned int i = 0;i < m_ComputationMask->GetImageDimension();++i)
        tmpSize[i] = maxPos[i] - minPos[i] + 1;

    m_ComputationRegion.SetSize(tmpSize);
}

template< typename TInputImage, typename TOutputImage >
void
MaskedImageToImageFilter < TInputImage, TOutputImage >
::BeforeThreadedGenerateData()
{
    if (m_ComputationRegion.GetSize(0) == 0)
        this->InitializeComputationRegionFromMask();

    unsigned int maxNumSlices = m_ComputationRegion.GetSize()[0];
    m_ProcessedDimension = 0;
    for (unsigned int i = 1;i < InputImageRegionType::ImageDimension;++i)
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
MaskedImageToImageFilter < TInputImage, TOutputImage >
::GenerateData()
{
    this->AllocateOutputs();
    this->BeforeThreadedGenerateData();

    this->InitializeSplitParameters();

    ThreadStruct str;
    str.Filter = this;

    this->GetMultiThreader()->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    this->GetMultiThreader()->SetSingleMethod(this->ThreaderMultiSplitCallback, &str);

    this->GetMultiThreader()->SingleMethodExecute();

    this->DeleteProgressReporter();

    this->AfterThreadedGenerateData();
}

template< typename TInputImage, typename TOutputImage >
void
MaskedImageToImageFilter < TInputImage, TOutputImage >
::InitializeSplitParameters()
{
    m_HighestProcessedSlice = 0;
    if (m_VerboseProgression)
    {
        if (m_ProgressReport)
            delete m_ProgressReport;

        m_ProgressReport = new itk::ProgressReporter(this,0,m_ComputationRegion.GetSize()[m_ProcessedDimension]);
    }
}

template< typename TInputImage, typename TOutputImage >
void
MaskedImageToImageFilter < TInputImage, TOutputImage >
::DeleteProgressReporter()
{
    if (m_ProgressReport)
    {
        delete m_ProgressReport;
        m_ProgressReport = 0;
    }
}

template< typename TInputImage, typename TOutputImage >
ITK_THREAD_RETURN_TYPE
MaskedImageToImageFilter < TInputImage, TOutputImage >
::ThreaderMultiSplitCallback(void *arg)
{
    ThreadStruct *str;
    itk::ThreadIdType threadId;

    threadId = ( (itk::MultiThreader::WorkUnitInfo *)( arg ) )->WorkUnitID;
    str = (ThreadStruct *)( ( (itk::MultiThreader::WorkUnitInfo *)( arg ) )->UserData );

    Self *filterPtr = dynamic_cast <Self *> (str->Filter.GetPointer());
    filterPtr->ThreadProcessSlices(threadId);

    return ITK_THREAD_RETURN_VALUE;
}

template< typename TInputImage, typename TOutputImage >
void
MaskedImageToImageFilter < TInputImage, TOutputImage >
::ThreadProcessSlices(itk::ThreadIdType threadId)
{
    InputImageRegionType processedRegion = m_ComputationRegion;
    processedRegion.SetSize(m_ProcessedDimension,1);
    unsigned int highestToleratedSliceValue = m_ComputationRegion.GetSize()[m_ProcessedDimension] - 1;

    bool continueLoop = true;
    while (continueLoop)
    {
        m_LockHighestProcessedSlice.Lock();

        if (m_HighestProcessedSlice > highestToleratedSliceValue)
        {
            m_LockHighestProcessedSlice.Unlock();
            continueLoop = false;
            continue;
        }

        processedRegion.SetIndex(m_ProcessedDimension, m_ComputationRegion.GetIndex()[m_ProcessedDimension] + m_HighestProcessedSlice);
        m_HighestProcessedSlice++;

        m_LockHighestProcessedSlice.Unlock();

        this->ThreadedGenerateData(processedRegion,threadId);

        if (m_VerboseProgression)
        {
            m_LockHighestProcessedSlice.Lock();
            m_ProgressReport->CompletedPixel();
            m_LockHighestProcessedSlice.Unlock();
        }
    }
}

} // end namespace anima
