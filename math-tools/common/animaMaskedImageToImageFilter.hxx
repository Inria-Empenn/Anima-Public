#pragma once
#include "animaMaskedImageToImageFilter.h"

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

        m_ComputationMask->SetRegions (this->GetInput(0)->GetLargestPossibleRegion());
        m_ComputationMask->SetSpacing (this->GetInput(0)->GetSpacing());
        m_ComputationMask->SetOrigin (this->GetInput(0)->GetOrigin());
        m_ComputationMask->SetDirection (this->GetInput(0)->GetDirection());

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

    MaskImageType::IndexType minPos, maxPos;

    for (unsigned int i = 0;i < m_ComputationMask->GetImageDimension();++i)
    {
        minPos[i] = m_ComputationMask->GetLargestPossibleRegion().GetIndex()[i] + m_ComputationMask->GetLargestPossibleRegion().GetSize()[i];
        maxPos[i] = 0;
    }

    while (!maskItr.IsAtEnd())
    {
        if (maskItr.Get() != 0)
        {
            MaskImageType::IndexType tmpInd = maskItr.GetIndex();

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

    MaskImageType::SizeType tmpSize;
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

    m_RegionSplitter->InitializeFromMask(m_ComputationMask,m_ComputationRegion,this->GetNumberOfThreads());

    // Since image requested region may now be smaller than the image, fill the outputs with zeros
    typedef typename TOutputImage::PixelType OutputPixelType;

    for (unsigned int i = 0;i < this->GetNumberOfOutputs();++i)
    {
        OutputPixelType zeroPixel;
        this->InitializeZeroPixel(this->GetOutput(i),zeroPixel);
        this->GetOutput(i)->FillBuffer(zeroPixel);
    }
}

} // end namespace anima
