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
::BeforeThreadedGenerateData()
{
    if (this->GetComputationRegion().GetSize(0) == 0)
        this->InitializeComputationRegionFromMask();

    Superclass::BeforeThreadedGenerateData();
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

    unsigned int numPoints = 0;
    while (!maskItr.IsAtEnd())
    {
        if (maskItr.Get() != 0)
        {
            ++numPoints;
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

    this->SetNumberOfPointsToProcess(numPoints);

    OutputImageRegionType computationRegion;
    computationRegion.SetIndex(minPos);

    MaskSizeType tmpSize;
    for (unsigned int i = 0;i < m_ComputationMask->GetImageDimension();++i)
        tmpSize[i] = maxPos[i] - minPos[i] + 1;

    computationRegion.SetSize(tmpSize);
    this->SetComputationRegion(computationRegion);
}

} // end namespace anima
