#pragma once
#include "animaStimulatedSpinEchoImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

#include <animaEPGSignalSimulator.h>

namespace anima
{
    
template <class TImage, class TOutputImage>
void
StimulatedSpinEchoImageFilter <TImage,TOutputImage>
::GenerateOutputInformation()
{
    this->Superclass::GenerateOutputInformation();

    OutputImageType *output = this->GetOutput();
    output->SetVectorLength(m_NumberOfEchoes);
}

template <class TImage, class TOutputImage>
StimulatedSpinEchoImageFilter <TImage,TOutputImage>
::StimulatedSpinEchoImageFilter()
{
    this->SetNumberOfRequiredInputs(3);

    m_NumberOfEchoes = 1;
    m_EchoSpacing = 10;
    m_ExcitationFlipAngle = M_PI / 2.0;
    m_FlipAngle = M_PI;
}

template <class TImage, class TOutputImage>
void StimulatedSpinEchoImageFilter <TImage,TOutputImage>
::SetInputT1(const TImage* T1)
{
    SetInput(0, const_cast<TImage*>(T1));
}

template <class TImage, class TOutputImage>
void StimulatedSpinEchoImageFilter <TImage,TOutputImage>
::SetInputT2(const TImage* T2)
{
    SetInput(1, const_cast<TImage*>(T2));
}

template <class TImage, class TOutputImage>
void StimulatedSpinEchoImageFilter <TImage,TOutputImage>
::SetInputM0(const TImage* M0)
{
    SetInput(2, const_cast<TImage*>(M0));
}

template <class TImage, class TOutputImage>
void StimulatedSpinEchoImageFilter <TImage,TOutputImage>
::SetInputB1(const TImage* B1)
{
    this->SetInput(3, const_cast<TImage*>(B1));
}

template <class TImage, class TOutputImage>
typename StimulatedSpinEchoImageFilter <TImage,TOutputImage>::Image4DType *
StimulatedSpinEchoImageFilter <TImage,TOutputImage>
::GetOutputAs4DImage()
{
    OutputImageType* outPtr = this->GetOutput();
    unsigned int ndim = outPtr->GetNumberOfComponentsPerPixel();

    typename Image4DType::RegionType region4d;
    typename TImage::RegionType region3d = outPtr->GetLargestPossibleRegion();
    typename Image4DType::SizeType size4d;
    typename Image4DType::SpacingType spacing4d;
    typename Image4DType::PointType origin4d;
    typename Image4DType::DirectionType direction4d;

    region4d.SetIndex(3,0);
    region4d.SetSize(3,ndim);
    size4d[3] = ndim;
    spacing4d[3] = 1;
    origin4d[3] = 0;
    direction4d(3,3) = 1;

    for (unsigned int i = 0;i < 3;++i)
    {
        region4d.SetIndex(i,region3d.GetIndex()[i]);
        region4d.SetSize(i,region3d.GetSize()[i]);

        size4d[i] = region3d.GetSize()[i];
        spacing4d[i] = outPtr->GetSpacing()[i];
        origin4d[i] = outPtr->GetOrigin()[i];
        for (unsigned int j = 0;j < 3;++j)
            direction4d(i,j) = (outPtr->GetDirection())(i,j);

        direction4d(i,3) = 0;
        direction4d(3,i) = 0;
    }

    m_Output4D = Image4DType::New();
    m_Output4D->Initialize();
    m_Output4D->SetRegions(region4d);
    m_Output4D->SetSpacing (spacing4d);
    m_Output4D->SetOrigin (origin4d);
    m_Output4D->SetDirection (direction4d);
    m_Output4D->Allocate();

    typedef itk::ImageRegionConstIterator <OutputImageType> OutImageIteratorType;
    typedef itk::ImageRegionIterator <Image4DType> Image4DIteratorType;

    for (unsigned int i = 0;i < ndim;++i)
    {
        region4d.SetIndex(3,i);
        region4d.SetSize(3,1);

        OutImageIteratorType outItr(outPtr,region3d);
        Image4DIteratorType fillItr(m_Output4D,region4d);

        while (!fillItr.IsAtEnd())
        {
            fillItr.Set(outItr.Get()[i]);
            ++fillItr;
            ++outItr;
        }
    }

    return m_Output4D;
}

template <class TImage, class TOutputImage>
void StimulatedSpinEchoImageFilter <TImage,TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread,
                       itk::ThreadIdType threadId)
{
    typename TImage::ConstPointer T1 = this->GetInput(0);
    typename TImage::ConstPointer T2 = this->GetInput(1);
    typename TImage::ConstPointer M0 = this->GetInput(2);
    typename TImage::ConstPointer B1;

    itk::ImageRegionIterator <TOutputImage> outputIterator(this->GetOutput(), outputRegionForThread);
    itk::ImageRegionConstIterator <TImage> inputIteratorT1(T1, outputRegionForThread);
    itk::ImageRegionConstIterator <TImage> inputIteratorT2(T2, outputRegionForThread);
    itk::ImageRegionConstIterator <TImage> inputIteratorM0(M0, outputRegionForThread);

    itk::ImageRegionConstIterator <TImage> inputIteratorB1;
    bool b1DataPresent = (this->GetNumberOfIndexedInputs() == 4);
    if (b1DataPresent)
    {
        B1 = this->GetInput(3);
        inputIteratorB1 = itk::ImageRegionConstIterator <TImage> (B1, outputRegionForThread);
    }

    typedef typename TOutputImage::PixelType VectorType;
    VectorType outputVector(m_NumberOfEchoes);

    anima::EPGSignalSimulator t2SignalSimulator;
    t2SignalSimulator.SetNumberOfEchoes(m_NumberOfEchoes);
    t2SignalSimulator.SetEchoSpacing(m_EchoSpacing);
    t2SignalSimulator.SetExcitationFlipAngle(m_ExcitationFlipAngle);

    anima::EPGSignalSimulator::RealVectorType tmpOutputVector;

    // For each voxel, calculate the signal of SP-GRE
    while(!outputIterator.IsAtEnd())
    {
        outputVector.Fill(0);
        if ((inputIteratorT1.Get() <= 0) || ( inputIteratorT2.Get() <= 0))
        {
            outputIterator.Set(outputVector);

            ++inputIteratorT1;
            ++inputIteratorT2;
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
        double t2Value = inputIteratorT2.Get();
        double m0Value = inputIteratorM0.Get();

        tmpOutputVector = t2SignalSimulator.GetValue(t1Value,t2Value,b1Value * m_FlipAngle,m0Value);

        for (unsigned int i = 0;i < m_NumberOfEchoes;++i)
            outputVector[i] = tmpOutputVector[i];

        outputIterator.Set(outputVector);

        ++inputIteratorT1;
        ++inputIteratorT2;
        ++inputIteratorM0;
        if (b1DataPresent)
            ++inputIteratorB1;

        ++outputIterator;
    }
}
    
}// end of namespace anima
