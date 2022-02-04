#pragma once
#include "animaCBFEstimationImageFilter_PCASL.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>

namespace anima
{

template <class InputPixelType, class OutputPixelType>
void
CBFEstimationImageFilter_PCASL <InputPixelType,OutputPixelType>
::BeforeThreadedGenerateData ()
{
    this->Superclass::BeforeThreadedGenerateData();

    if (!m_M0Image)
    {
        m_M0Image = InputImageType::New();
        m_M0Image->Initialize();
        m_M0Image->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
        m_M0Image->SetSpacing (this->GetInput(0)->GetSpacing());
        m_M0Image->SetOrigin (this->GetInput(0)->GetOrigin());
        m_M0Image->SetDirection (this->GetInput(0)->GetDirection());
        m_M0Image->Allocate();
        m_M0Image->FillBuffer(m_M0ConstantValue);
    }
}

template <class InputPixelType, class OutputPixelType>
void
CBFEstimationImageFilter_PCASL <InputPixelType,OutputPixelType>
::ThreadedGenerateData (const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionConstIteratorWithIndex <InputImageType> InputImageIteratorType;
    typedef itk::ImageRegionConstIterator <MaskImageType> MaskIteratorType;
    typedef itk::ImageRegionIterator <OutputImageType> OutputImageIteratorType;
    
    InputImageIteratorType inputItr(this->GetInput(0), outputRegionForThread);
    InputImageIteratorType m0Itr(m_M0Image, outputRegionForThread);
    MaskIteratorType maskItr(this->GetComputationMask(),outputRegionForThread);
    OutputImageIteratorType outItr(this->GetOutput(0),outputRegionForThread);
    
    IndexType currentIndex;
    double logConstantValue = std::log(6.0e6) + std::log(m_LambdaParameter) - std::log(m_AlphaParameter) - std::log(m_BloodT1) - std::log(1.0 - std::exp(- m_LabelDuration / m_BloodT1));

    while (!inputItr.IsAtEnd())
    {
        if (maskItr.Get() == 0)
        {
            ++maskItr;
            ++m0Itr;
            ++inputItr;
            ++outItr;

            continue;
        }

        currentIndex = inputItr.GetIndex();
        double currentPLD = currentIndex[InputImageType::ImageDimension - 1] * m_SliceDelay + m_BasePostLabelingDelay;

        double logM0 = std::log(1.0e-16);
        if (m0Itr.Get() > 1.0e-16)
            logM0 = std::log(m0Itr.Get());

        double logInput = std::log(1.0e-16);
        if (inputItr.Get() > 1.0e-16)
            logInput = std::log(inputItr.Get());

        double logCBFValue = logConstantValue + logInput - logM0 + currentPLD / m_BloodT1;

        double cbfValue = std::exp(logCBFValue);
        if (!std::isfinite(cbfValue))
            cbfValue = 0;

        outItr.Set(cbfValue);
        
        ++maskItr;
        ++m0Itr;
        ++inputItr;
        ++outItr;
    }
}

} //end of namespace anima

