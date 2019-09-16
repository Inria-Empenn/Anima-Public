#pragma once
#include "animaSimuBlochSP-GRE.h"

#include <itkObjectFactory.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

namespace anima
{

template< class TImage>
SimuBlochSPGRE<TImage>::SimuBlochSPGRE()
{
    this->SetNumberOfRequiredInputs(3);
}

template< class TImage>
void SimuBlochSPGRE<TImage>::SetInputT1(const TImage* T1)
{
    SetInput(0, const_cast<TImage*>(T1));
}

template< class TImage>
void SimuBlochSPGRE<TImage>::SetInputT2s(const TImage* T2s)
{
    SetInput(1, const_cast<TImage*>(T2s));
}

template< class TImage>
void SimuBlochSPGRE<TImage>::SetInputM0(const TImage* M0)
{
    SetInput(2, const_cast<TImage*>(M0));
}

template< class TImage>
void SimuBlochSPGRE<TImage>::SetInputB1(const TImage* B1)
{
    this->SetInput(3, const_cast<TImage*>(B1));
}

template< class TImage>
void SimuBlochSPGRE< TImage>
::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    typename TImage::ConstPointer T1 = this->GetInput(0);
    typename TImage::ConstPointer T2s = this->GetInput(1);
    typename TImage::ConstPointer M0 = this->GetInput(2);
    typename TImage::ConstPointer B1;

    itk::ImageRegionIterator<TImage> outputIterator(this->GetOutput(), outputRegionForThread);
    itk::ImageRegionConstIterator<TImage> inputIteratorT1(T1, outputRegionForThread);
    itk::ImageRegionConstIterator<TImage> inputIteratorT2s(T2s, outputRegionForThread);
    itk::ImageRegionConstIterator<TImage> inputIteratorM0(M0, outputRegionForThread);

    itk::ImageRegionConstIterator<TImage> inputIteratorB1;
    bool b1DataPresent = (this->GetNumberOfIndexedInputs() == 4);
    if (b1DataPresent)
    {
        B1 = this->GetInput(3);
        inputIteratorB1 = itk::ImageRegionConstIterator<TImage> (B1, outputRegionForThread);
    }

    // For each voxel, calculate the signal of SP-GRE
    while(!outputIterator.IsAtEnd())
    {
        if ((inputIteratorT1.Get() <= 0) || (inputIteratorT2s.Get() <= 0))
        {
            outputIterator.Set(0);

            ++inputIteratorT1;
            ++inputIteratorT2s;
            ++inputIteratorM0;
            if (b1DataPresent)
                ++inputIteratorB1;

            ++outputIterator;
            continue;
        }

        double b1Value = 1;
        if (b1DataPresent)
            b1Value = inputIteratorB1.Get();

        double t1Value = inputIteratorT1.Get();
        double t2sValue = inputIteratorT2s.Get();
        double m0Value = inputIteratorM0.Get();

        double sinFA = std::sin(M_PI * b1Value * m_FA / 180);
        double cosFA = std::cos(M_PI * b1Value * m_FA / 180);

        outputIterator.Set(m0Value * (1- exp(- m_TR / t1Value)) * sinFA / (1 - exp(- m_TR / t1Value) * cosFA ) * exp(- m_TE / t2sValue));

        ++inputIteratorT1;
        ++inputIteratorT2s;
        ++inputIteratorM0;
        if (b1DataPresent)
            ++inputIteratorB1;

        ++outputIterator;
    }
}
    
}// end of namespace anima
