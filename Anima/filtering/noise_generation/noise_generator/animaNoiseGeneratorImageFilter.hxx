#pragma once

#include "animaNoiseGeneratorImageFilter.h"
#include <animaDistributionSampling.h>
#include <animaKummerFunctions.h>

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

namespace anima
{

template <unsigned int Dimension>
void
NoiseGeneratorImageFilter<Dimension>
::GenerateOutputInformation ()
{
    //-------------------------------
    // Creates the additional outputs
    //-------------------------------
    unsigned int prevNum = this->GetNumberOfOutputs();
    
    this->SetNumberOfIndexedOutputs(m_NumberOfReplicates);
    
    for (unsigned int i = prevNum;i < m_NumberOfReplicates;++i)
    {
        this->SetNthOutput(i,this->MakeOutput(i).GetPointer());
    }
    
    //---------------------------------------------------
    // Call the superclass' implementation of this method
    //---------------------------------------------------
    Superclass::GenerateOutputInformation();
}

template <unsigned int Dimension>
void
NoiseGeneratorImageFilter<Dimension>
::BeforeThreadedGenerateData ()
{
    m_Generators.clear();
    std::mt19937 motherGenerator(time(0));
    for (int i = 0;i < this->GetNumberOfThreads();++i)
        m_Generators.push_back(std::mt19937(motherGenerator()));
}

template <unsigned int Dimension>
void
NoiseGeneratorImageFilter<Dimension>
::ThreadedGenerateData (const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    if (Dimension < 4)
        this->TreatRegionWithNoiseVariance(outputRegionForThread,threadId);
    else
    {
        OutputImageRegionType region = outputRegionForThread;
        region.SetSize(3,1);
        unsigned int maxIndex4d = outputRegionForThread.GetIndex()[3] + outputRegionForThread.GetSize()[3];

        for (unsigned int i = outputRegionForThread.GetIndex()[3];i < maxIndex4d;++i)
        {
            region.SetIndex(3,i);
            this->TreatRegionWithNoiseVariance(region,threadId);
        }
    }
}

template <unsigned int Dimension>
void
NoiseGeneratorImageFilter<Dimension>
::TreatRegionWithNoiseVariance(const OutputImageRegionType &region, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionConstIterator<TInputImage> InputImageIteratorType;
    typedef itk::ImageRegionIterator<OutputImageType> OutputImageIteratorType;

    InputImageIteratorType inputIterator(this->GetInput(), region);
    
    std::vector<OutputImageIteratorType> outIterators(m_NumberOfReplicates);
    for (unsigned int i = 0;i < m_NumberOfReplicates;++i)
        outIterators[i] = OutputImageIteratorType(this->GetOutput(i), region);
    
    double mean = 0;

    while (!inputIterator.IsAtEnd())
    {
        double refData = inputIterator.Get();
        double sd = m_StandardDeviation;
        
        if (m_UseGaussianDistribution)
        {
            double varianceValue = m_StandardDeviation * m_StandardDeviation;
            double sqSignal = refData * refData;
            double x = -1.0 * sqSignal / (2.0 * varianceValue);
            double laguerreValue = KummerFunction(x, -0.5, 1.0);
            sd = std::sqrt( 2.0 * varianceValue + sqSignal - M_PI * varianceValue * laguerreValue * laguerreValue / 2.0 );
            
            std::cout << m_StandardDeviation << " " << sd << std::endl;
        }
        
        for (unsigned int i = 0;i < m_NumberOfReplicates;++i)
        {
            double realNoise = SampleFromGaussianDistribution(mean, sd, m_Generators[threadId]);
            double data = refData + realNoise;
            
            if (!m_UseGaussianDistribution)
            {
                double imagNoise = SampleFromGaussianDistribution(mean, sd, m_Generators[threadId]);
                data = std::sqrt(data * data + imagNoise * imagNoise);
            }
            
            if ((std::isnan(data)) || (!std::isfinite(data)))
                data = refData;
            
            outIterators[i].Set(data);
            ++outIterators[i];
        }

        ++inputIterator;
    }
}

} //end of namespace anima
