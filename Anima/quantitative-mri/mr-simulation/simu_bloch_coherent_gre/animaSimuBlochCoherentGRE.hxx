#pragma once
#include "animaSimuBlochCoherentGRE.h"

#include <itkObjectFactory.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

namespace anima
{

template< class TImage>
SimuBlochCoherentGRE<TImage>::SimuBlochCoherentGRE()
{
    this->SetNumberOfRequiredInputs(3);
}

template< class TImage>
void SimuBlochCoherentGRE<TImage>::SetInputT1(const TImage* T1)
{
    SetInput(0, const_cast<TImage*>(T1));
}

template< class TImage>
void SimuBlochCoherentGRE<TImage>::SetInputT2s(const TImage* T2s)
{
    SetInput(1, const_cast<TImage*>(T2s));
}

template< class TImage>
void SimuBlochCoherentGRE<TImage>::SetInputM0(const TImage* M0)
{
    SetInput(2, const_cast<TImage*>(M0));
}

template< class TImage>
void SimuBlochCoherentGRE<TImage>::SetInputT2(const TImage* T2)//changed for CoherentGRE
{
    SetInput(3, const_cast<TImage*>(T2));//changed for CoherentGRE
}

template< class TImage>//changed for CoherentGRE

void SimuBlochCoherentGRE< TImage>
::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    typename TImage::ConstPointer T1 = this->GetInput(0);
    typename TImage::ConstPointer T2s = this->GetInput(1);
    typename TImage::ConstPointer M0 = this->GetInput(2);
    typename TImage::ConstPointer T2 = this->GetInput(3);//changed for CoherentGRE

    itk::ImageRegionIterator<TImage> outputIterator(this->GetOutput(), outputRegionForThread);
    itk::ImageRegionConstIterator<TImage> inputIteratorT1(T1, outputRegionForThread);
    itk::ImageRegionConstIterator<TImage> inputIteratorT2s(T2s, outputRegionForThread);
    itk::ImageRegionConstIterator<TImage> inputIteratorM0(M0, outputRegionForThread);
    itk::ImageRegionConstIterator<TImage> inputIteratorT2(T2, outputRegionForThread);//changed for CoherentGRE

    double sinFA= sin(M_PI * m_FA / 180);
    double cosFA = cos(M_PI * m_FA / 180);

    // For each voxel, calculate the signal of CoherentGRE
    while(!outputIterator.IsAtEnd())
    {
        if ( (inputIteratorT1.Get() <= 0) || ( inputIteratorT2s.Get() <=0) || ( inputIteratorT2.Get() <=0))//changed for CoherentGRE
        {
            outputIterator.Set(0);
        }
        else
        {
            //      outputIterator.Set(inputIteratorM0.Get() * (1- exp(- m_TR / inputIteratorT1.Get())) * sinFA / (1 - exp(- m_TR / inputIteratorT1.Get()) * cosFA  - exp( - m_TR /  inputIteratorT2.Get() ) * (exp( - m_TR / inputIteratorT1.Get() ) - cosFA) ) * exp(- m_TE / inputIteratorT2s.Get()) ); //changed for CoherentGRE//changed for CoherentGRE
            outputIterator.Set(inputIteratorM0.Get() * sinFA / (1 + inputIteratorT1.Get() / inputIteratorT2.Get() - cosFA * ( inputIteratorT1.Get()/ inputIteratorT2.Get() -1 ) ) * exp(- m_TE / inputIteratorT2s.Get()) ); //changed for CoherentGRE//changed for CoherentGRE
        }

        ++inputIteratorT1;
        ++inputIteratorT2s;
        ++inputIteratorM0;
        ++inputIteratorT2;//changed for CoherentGRE

        ++outputIterator;
    }
}
    
}// end of namespace anima
