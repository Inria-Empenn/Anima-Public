#pragma once
#include "animaSimuBlochIR-SE.h"

#include <itkObjectFactory.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

namespace anima
{

template< class TImage>
SimuBlochIRSE<TImage>::SimuBlochIRSE()
{
    this->SetNumberOfRequiredInputs(3);
}

template< class TImage>
void SimuBlochIRSE<TImage>::SetInputT1(const TImage* T1)
{
    SetInput(0, const_cast<TImage*>(T1));
}

template< class TImage>
void SimuBlochIRSE<TImage>::SetInputT2(const TImage* T2)
{
    SetInput(1, const_cast<TImage*>(T2));
}

template< class TImage>
void SimuBlochIRSE<TImage>::SetInputM0(const TImage* M0)
{
    SetInput(2, const_cast<TImage*>(M0));
}

template< class TImage>
void SimuBlochIRSE< TImage>
::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    typename TImage::ConstPointer T1 = this->GetInput(0);
    typename TImage::ConstPointer T2 = this->GetInput(1);
    typename TImage::ConstPointer M0 = this->GetInput(2);

    itk::ImageRegionIterator<TImage> outputIterator(this->GetOutput(), outputRegionForThread);
    itk::ImageRegionConstIterator<TImage> inputIteratorT1(T1, outputRegionForThread);
    itk::ImageRegionConstIterator<TImage> inputIteratorT2(T2, outputRegionForThread);
    itk::ImageRegionConstIterator<TImage> inputIteratorM0(M0, outputRegionForThread);

    // For each voxel, calculate the signal of Inversion Recovery using bloch equation
    while(!outputIterator.IsAtEnd())
    {
        if ( (inputIteratorT1.Get() <= 0) || ( inputIteratorT2.Get() <= 0))
        {
            outputIterator.Set(0);
        }
        else
        {
            outputIterator.Set(inputIteratorM0.Get() * fabs(1 - 2 * exp(- m_TI / inputIteratorT1.Get()) + exp(- m_TR / inputIteratorT1.Get())) * exp(- m_TE / inputIteratorT2.Get()) );
        }

        ++inputIteratorT1;
        ++inputIteratorT2;
        ++inputIteratorM0;

        ++outputIterator;
    }
}
    
}// end of namespace anima
