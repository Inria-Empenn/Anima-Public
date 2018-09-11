#pragma once
#include "animaMaskedImageToImageFilter.h"

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
    this->CheckComputationMask();
    Superclass::BeforeThreadedGenerateData();

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
