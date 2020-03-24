#pragma once

#include "animaNoiseGeneratorImageFilter.h"
#include <animaDistributionSampling.h>

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
    Superclass::BeforeThreadedGenerateData();

    m_Generators.clear();
    std::mt19937 motherGenerator(time(ITK_NULLPTR));
    for (int i = 0;i < this->GetNumberOfWorkUnits();++i)
        m_Generators.push_back(std::mt19937(motherGenerator()));
}

template <class ImageType>
void
NoiseGeneratorImageFilter<ImageType>
::DynamicThreadedGenerateData (const OutputImageRegionType &outputRegionForThread)
{
    typedef itk::ImageRegionConstIterator<InputImageType> InputImageIteratorType;
    typedef itk::ImageRegionIterator<OutputImageType> OutputImageIteratorType;
    
    InputImageIteratorType inputIterator(this->GetInput(), outputRegionForThread);
    
    std::vector<OutputImageIteratorType> outIterators(m_NumberOfReplicates);
    for (unsigned int i = 0;i < m_NumberOfReplicates;++i)
        outIterators[i] = OutputImageIteratorType(this->GetOutput(i), outputRegionForThread);
    
    unsigned int threadId = this->GetSafeThreadId();
    
    while (!inputIterator.IsAtEnd())
    {
        double refData = inputIterator.Get();
        
        for (unsigned int i = 0;i < m_NumberOfReplicates;++i)
        {
            double realNoise = SampleFromGaussianDistribution(0.0, m_NoiseSigma, m_Generators[threadId]);
            double data = refData + realNoise;
            
            if (!m_UseGaussianDistribution)
            {
                double imagNoise = SampleFromGaussianDistribution(0.0, m_NoiseSigma, m_Generators[threadId]);
                data = std::sqrt(data * data + imagNoise * imagNoise);
            }
            
            if ((std::isnan(data)) || (!std::isfinite(data)))
                data = refData;
            
            outIterators[i].Set(data);
            ++outIterators[i];
        }
        
        ++inputIterator;
    }

    this->SafeReleaseThreadId(threadId);
}

} //end of namespace anima
