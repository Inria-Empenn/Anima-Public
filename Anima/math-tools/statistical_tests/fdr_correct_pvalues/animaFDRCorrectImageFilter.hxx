#pragma once
#include "animaFDRCorrectImageFilter.h"

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

#include <animaFDRCorrection.h>

namespace anima
{

template <class PixelScalarType>
void
FDRCorrectImageFilter <PixelScalarType>
::GenerateData()
{
    this->AllocateOutputs();

    if (!m_MaskImage)
    {
        this->CreateFullMask();
        this->GetOutput()->SetRequestedRegion(this->GetInput()->GetLargestPossibleRegion());
    }

    typedef itk::ImageRegionConstIterator <TInputImage> InputImageIteratorType;
    typedef itk::ImageRegionConstIterator <MaskImageType> MaskIteratorType;
    typedef itk::ImageRegionIterator <TOutputImage> OutputImageIteratorType;

    MaskIteratorType maskItr(m_MaskImage, this->GetOutput()->GetRequestedRegion());
    InputImageIteratorType inItr(this->GetInput(), this->GetOutput()->GetRequestedRegion());
    OutputImageIteratorType outItr(this->GetOutput(), this->GetOutput()->GetRequestedRegion());

    std::vector <double> pvalues;

    while (!maskItr.IsAtEnd())
    {
        if (maskItr.Get() == 0)
        {
            ++maskItr;
            ++inItr;
            continue;
        }

        pvalues.push_back(inItr.Get());

        ++inItr;
        ++maskItr;
    }

    anima::BYCorrection(pvalues, m_QValue);

    maskItr.GoToBegin();

    unsigned int pos = 0;
    while (!maskItr.IsAtEnd())
    {
        if (maskItr.Get() == 0)
        {
            ++maskItr;
            outItr.Set(0);
            ++outItr;
            continue;
        }

        outItr.Set(pvalues[pos]);
        ++pos;

        ++outItr;
        ++maskItr;
    }
}

template <class PixelScalarType>
void
FDRCorrectImageFilter <PixelScalarType>
::CreateFullMask()
{
    if (m_MaskImage)
        return;

    m_MaskImage = MaskImageType::New();
    m_MaskImage->Initialize();
    m_MaskImage->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    m_MaskImage->SetOrigin(this->GetInput()->GetOrigin());
    m_MaskImage->SetDirection(this->GetInput()->GetDirection());
    m_MaskImage->SetSpacing(this->GetInput()->GetSpacing());

    m_MaskImage->Allocate();

    m_MaskImage->FillBuffer(1);
}

} // end namespace anima
