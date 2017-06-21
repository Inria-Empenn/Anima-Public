#pragma once
#include "animaSVFLieBracketImageFilter.h"

#include <animaJacobianMatrixImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

namespace anima
{

template <typename TPixelType, unsigned int Dimension>
void
SVFLieBracketImageFilter <TPixelType, Dimension>
::BeforeThreadedGenerateData()
{
    this->Superclass::BeforeThreadedGenerateData();

    typedef anima::JacobianMatrixImageFilter <TPixelType, TPixelType, Dimension> JacobianFilterType;

    if (m_FirstFieldJacobian.IsNull())
    {
        typename JacobianFilterType::Pointer jacFilter = JacobianFilterType::New();
        jacFilter->SetInput(this->GetInput(0));
        jacFilter->SetNoIdentity(true);
        jacFilter->SetNumberOfThreads(this->GetNumberOfThreads());
        jacFilter->SetNeighborhood(1);

        jacFilter->Update();

        m_FirstFieldJacobian = jacFilter->GetOutput();
        m_FirstFieldJacobian->DisconnectPipeline();
    }

    if (m_SecondFieldJacobian.IsNull())
    {
        typename JacobianFilterType::Pointer jacFilter = JacobianFilterType::New();
        jacFilter->SetInput(this->GetInput(1));
        jacFilter->SetNoIdentity(true);
        jacFilter->SetNumberOfThreads(this->GetNumberOfThreads());
        jacFilter->SetNeighborhood(1);

        jacFilter->Update();

        m_SecondFieldJacobian = jacFilter->GetOutput();
        m_SecondFieldJacobian->DisconnectPipeline();
    }
}

template <typename TPixelType, unsigned int Dimension>
void
SVFLieBracketImageFilter <TPixelType, Dimension>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionConstIterator <InputImageType> InputIteratorType;
    typedef itk::ImageRegionConstIterator <JacobianImageType> JacobianIteratorType;
    typedef itk::ImageRegionIterator <OutputImageType> OutIteratorType;

    InputIteratorType firstInputItr(this->GetInput(0),outputRegionForThread);
    InputIteratorType secondInputItr(this->GetInput(1),outputRegionForThread);
    OutIteratorType outItr(this->GetOutput(),outputRegionForThread);

    JacobianIteratorType firstJacItr(m_FirstFieldJacobian,outputRegionForThread);
    JacobianIteratorType secondJacItr(m_SecondFieldJacobian,outputRegionForThread);

    InputPixelType inputValue;
    OutputPixelType outputValue;
    JacobianPixelType jacValue;

    while (!outItr.IsAtEnd())
    {
        outputValue.Fill(0);

        inputValue = secondInputItr.Get();
        jacValue = firstJacItr.Get();

        for (unsigned int i = 0;i < Dimension;++i)
        {
            for (unsigned int j = 0;j < Dimension;++j)
                outputValue[i] += jacValue[i * Dimension + j] * inputValue[j];
        }

        inputValue = firstInputItr.Get();
        jacValue = secondJacItr.Get();

        for (unsigned int i = 0;i < Dimension;++i)
        {
            for (unsigned int j = 0;j < Dimension;++j)
                outputValue[i] -= jacValue[i * Dimension + j] * inputValue[j];
        }

        outItr.Set(outputValue);

        ++firstInputItr;
        ++secondInputItr;
        ++outItr;
        ++firstJacItr;
        ++secondJacItr;
    }
}

} // end namespace anima
