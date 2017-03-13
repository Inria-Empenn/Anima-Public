#pragma once
#include "animaT1RelaxometryEstimationImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

namespace anima
{
template <typename TInputImage, typename TOutputImage>
void
T1RelaxometryEstimationImageFilter <TInputImage,TOutputImage>
::BeforeThreadedGenerateData()
{
    if (this->GetNumberOfIndexedInputs() != m_FlipAngles.size())
        itkExceptionMacro("There should be the same number of inputs and flip angles");

    if (this->GetNumberOfIndexedInputs() != 2)
        itkExceptionMacro("There should be only two inputs for DESPOT1");

    Superclass::BeforeThreadedGenerateData();
}

template <typename TInputImage, typename TOutputImage>
void
T1RelaxometryEstimationImageFilter <TInputImage,TOutputImage>
::CheckComputationMask()
{
    if (this->GetComputationMask())
        return;

    typedef itk::ImageRegionConstIterator <InputImageType> IteratorType;
    typedef itk::ImageRegionIterator <MaskImageType> MaskIteratorType;

    std::vector <IteratorType> inItrs(this->GetNumberOfIndexedInputs());
    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
        inItrs[i] = IteratorType(this->GetInput(i),this->GetOutput()->GetLargestPossibleRegion());

    typename MaskImageType::Pointer maskImage = MaskImageType::New();
    maskImage->Initialize();
    maskImage->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
    maskImage->SetSpacing (this->GetInput(0)->GetSpacing());
    maskImage->SetOrigin (this->GetInput(0)->GetOrigin());
    maskImage->SetDirection (this->GetInput(0)->GetDirection());
    maskImage->Allocate();

    MaskIteratorType maskItr (maskImage,this->GetOutput()->GetLargestPossibleRegion());
    while (!maskItr.IsAtEnd())
    {
        double averageVal = 0;
        for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
            averageVal += inItrs[i].Get();

        averageVal /= this->GetNumberOfIndexedInputs();

        bool maskPoint = (averageVal <= m_AverageSignalThreshold);

        if (maskPoint)
            maskItr.Set(0);
        else
            maskItr.Set(1);

        for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
            ++inItrs[i];

        ++maskItr;
    }

    this->SetComputationMask(maskImage);
}

template <typename TInputImage, typename TOutputImage>
void
T1RelaxometryEstimationImageFilter <TInputImage,TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionConstIterator <InputImageType> ImageIteratorType;
    typedef itk::ImageRegionIterator <OutputImageType> OutImageIteratorType;

    OutImageIteratorType outT1Iterator(this->GetOutput(0),outputRegionForThread);
    OutImageIteratorType outM0Iterator(this->GetOutput(1),outputRegionForThread);

    ImageIteratorType firstInputIterator(this->GetInput(0),outputRegionForThread);
    ImageIteratorType secondInputIterator(this->GetInput(1),outputRegionForThread);

    OutImageIteratorType b1MapItr;
    if (m_B1Map)
        b1MapItr = OutImageIteratorType (m_B1Map,outputRegionForThread);

    typedef itk::ImageRegionIterator <MaskImageType> MaskIteratorType;
    MaskIteratorType maskItr(this->GetComputationMask(),outputRegionForThread);

    while (!maskItr.IsAtEnd())
    {
        double b1Value = 1;
        if (m_B1Map)
            b1Value = b1MapItr.Get();

        if ((maskItr.Get() == 0)||(b1Value == 0))
        {
            outT1Iterator.Set(0);
            outM0Iterator.Set(0);

            ++maskItr;
            ++outT1Iterator;
            ++outM0Iterator;

            ++firstInputIterator;
            ++secondInputIterator;

            if (m_B1Map)
                ++b1MapItr;

            continue;
        }

        double basicT1Value = m_TRValue;

        double lnValue = secondInputIterator.Get() * std::sin(m_FlipAngles[0] * b1Value) * std::cos(m_FlipAngles[1] * b1Value) - firstInputIterator.Get() * std::cos(m_FlipAngles[0] * b1Value) * std::sin(m_FlipAngles[1] * b1Value);
        lnValue /= secondInputIterator.Get() * std::sin(m_FlipAngles[0] * b1Value) - firstInputIterator.Get() * std::sin(m_FlipAngles[1] * b1Value);

        if (lnValue != lnValue) // test NaN
        {
            outT1Iterator.Set(1e-4);
            outM0Iterator.Set(0);

            ++maskItr;
            ++outT1Iterator;
            ++outM0Iterator;

            ++firstInputIterator;
            ++secondInputIterator;

            if (m_B1Map)
                ++b1MapItr;

            continue;
        }
        else
            lnValue = std::log(std::abs(lnValue));

        basicT1Value /= lnValue;

        if ((basicT1Value < 0)||(basicT1Value > m_T1UpperBoundValue))
            basicT1Value = m_T1UpperBoundValue;

        outT1Iterator.Set(basicT1Value);

        // M0
        double sumSignalData = 0;
        double sumSquaredData = 0;

        double dataValue = (1.0 - std::exp(- m_TRValue / basicT1Value)) * std::sin(m_FlipAngles[0] * b1Value) / (1.0 - std::exp(-m_TRValue / basicT1Value) * std::cos(m_FlipAngles[0] * b1Value));

        sumSignalData += firstInputIterator.Get() * dataValue;
        sumSquaredData += dataValue * dataValue;

        dataValue = (1.0 - std::exp(- m_TRValue / basicT1Value)) * std::sin(m_FlipAngles[1] * b1Value) / (1.0 - std::exp(- m_TRValue / basicT1Value) * std::cos(m_FlipAngles[1] * b1Value));

        sumSignalData += secondInputIterator.Get() * dataValue;
        sumSquaredData += dataValue * dataValue;

        double m0Value = sumSignalData / sumSquaredData;

        if (m0Value <= 1.0e-4)
            m0Value = 1.0e-4;
        else if (m0Value > m_M0UpperBoundValue)
            m0Value = m_M0UpperBoundValue;

        outM0Iterator.Set(m0Value);

        ++maskItr;
        ++outT1Iterator;
        ++outM0Iterator;

        ++firstInputIterator;
        ++secondInputIterator;

        if (m_B1Map)
            ++b1MapItr;
    }
}

} // end namespace anima
