#pragma once

#include "animaNoiseGeneratorImageFilter.h"
#include <animaDistributionSampling.h>
#include <animaKummerFunctions.h>

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

namespace anima
{

template <class ImageType>
void
NoiseGeneratorImageFilter<ImageType>
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

template <class ImageType>
void
NoiseGeneratorImageFilter<ImageType>
::BeforeThreadedGenerateData ()
{
    m_Generators.clear();
    std::mt19937 motherGenerator(time(0));
    for (int i = 0;i < this->GetNumberOfThreads();++i)
        m_Generators.push_back(std::mt19937(motherGenerator()));
}

template <class ImageType>
void
NoiseGeneratorImageFilter<ImageType>
::ThreadedGenerateData (const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionConstIterator<InputImageType> InputImageIteratorType;
    typedef itk::ImageRegionIterator<OutputImageType> OutputImageIteratorType;
    
    InputImageIteratorType inputIterator(this->GetInput(), outputRegionForThread);
    
    std::vector<OutputImageIteratorType> outIterators(m_NumberOfReplicates);
    for (unsigned int i = 0;i < m_NumberOfReplicates;++i)
        outIterators[i] = OutputImageIteratorType(this->GetOutput(i), outputRegionForThread);
    
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
