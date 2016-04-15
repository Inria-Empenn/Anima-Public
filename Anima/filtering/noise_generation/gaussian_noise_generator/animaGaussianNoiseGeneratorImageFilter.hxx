#pragma once

#include "animaGaussianNoiseGeneratorImageFilter.h"
#include <animaDistributionSampling.h>

#include <itkSymmetricEigenAnalysis.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>

namespace anima
{

template <unsigned int Dimension>
void
GaussianNoiseGeneratorImageFilter<Dimension>
::BeforeThreadedGenerateData ()
{
    unsigned int nbInputs = this->GetNumberOfIndexedInputs();
    if (nbInputs <= 0)
    {
        std::cerr << "Error: No inputs available... Exiting..." << std::endl;
        exit(-1);
    }

    m_MeanVals.clear();
    m_MeanVals.push_back(m_RefVal);

    if ((Dimension == 4)&&(!m_AverageMeanValOnAllVolumes))
    {
        for (unsigned int i = 1;i < this->GetInput()->GetLargestPossibleRegion().GetSize()[3];++i)
            m_MeanVals.push_back(m_RefVal);
    }

    if (m_RefVal == 0)
    {
        typedef itk::ImageRegionConstIteratorWithIndex< TInputImage > InIteratorType;

        if ((Dimension < 4)||(m_AverageMeanValOnAllVolumes))
        {
            InIteratorType inputIterator (this->GetInput(), this->GetInput()->GetLargestPossibleRegion());

            double meanVal = 0;
            unsigned int nbPts = 0;
            while (!inputIterator.IsAtEnd())
            {
                if (inputIterator.Get() < m_BackgroundThreshold)
                {
                    ++inputIterator;
                    continue;
                }

                meanVal += inputIterator.Get();
                ++nbPts;
                ++inputIterator;
            }

            m_MeanVals[0] = meanVal / nbPts;

            std::cout << "Reference image value: " << m_MeanVals[0] << ", noise stdev: "
            << m_RelStdDevGaussianNoise * m_MeanVals[0] << std::endl;
        }
        else
        {
            // Treat volumes independently
            OutputImageRegionType region = this->GetInput()->GetLargestPossibleRegion();
            region.SetSize(3,1);
            for (unsigned int i = 0;i < this->GetInput()->GetLargestPossibleRegion().GetSize()[3];++i)
            {
                region.SetIndex(3,i);
                InIteratorType inputIterator (this->GetInput(), region);

                double meanVal = 0;
                unsigned int nbPts = 0;
                while (!inputIterator.IsAtEnd())
                {
                    if (inputIterator.Get() < m_BackgroundThreshold)
                    {
                        ++inputIterator;
                        continue;
                    }

                    meanVal += inputIterator.Get();
                    ++nbPts;
                    ++inputIterator;
                }

                m_MeanVals[i] = meanVal / nbPts;
            }
        }
    }
    else
    {
        std::cout << "Reference image value: " << m_MeanVals[0] << ", noise stdev: "
        << m_RelStdDevGaussianNoise * m_MeanVals[0] << std::endl;
    }

    m_Generators.clear();

    std::mt19937 motherGenerator(time(0));

    for (int i = 0;i < this->GetNumberOfThreads();++i)
        m_Generators.push_back(std::mt19937(motherGenerator()));
}

template <unsigned int Dimension>
void
GaussianNoiseGeneratorImageFilter<Dimension>
::ThreadedGenerateData (const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    if ((Dimension < 4)||(m_AverageMeanValOnAllVolumes))
    {
        double varVal = m_RelStdDevGaussianNoise * m_MeanVals[0];

        this->TreatRegionWithNoiseVariance(outputRegionForThread,varVal,threadId);
    }
    else
    {
        OutputImageRegionType region = outputRegionForThread;
        region.SetSize(3,1);
        unsigned int maxIndex4d = outputRegionForThread.GetIndex()[3] + outputRegionForThread.GetSize()[3];

        for (unsigned int i = outputRegionForThread.GetIndex()[3];i < maxIndex4d;++i)
        {
            region.SetIndex(3,i);
            double varVal = m_RelStdDevGaussianNoise * m_MeanVals[i];

            this->TreatRegionWithNoiseVariance(region,varVal,threadId);
        }
    }
}

template <unsigned int Dimension>
void
GaussianNoiseGeneratorImageFilter<Dimension>
::TreatRegionWithNoiseVariance(const OutputImageRegionType &region, double &variance, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionConstIteratorWithIndex< TInputImage > InIteratorType;
    typedef itk::ImageRegionIteratorWithIndex< OutputImageType > OutRegionIteratorType;

    OutRegionIteratorType outIterator(this->GetOutput(), region);
    InIteratorType inputIterator (this->GetInput(), region);
    double mean = 0;

    while (!outIterator.IsAtEnd())
    {
        double refData = inputIterator.Get();

        double gaussNoise = anima::SampleFromGaussianDistribution(mean, variance, m_Generators[threadId]);
        double data = refData + gaussNoise;

        if ((std::isnan(data))||(!std::isfinite(data)))
            data = refData;

        outIterator.Set(data);

        ++inputIterator;
        ++outIterator;
    }
}

} //end of namespace anima
