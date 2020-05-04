#pragma once
#include "animaNonLocalMeansImageFilter.h"

#include <animaMeanAndVarianceImagesFilter.h>
#include <animaNonLocalMeansPatchSearcher.h>

#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>

namespace anima
{

template <class TInputImage>
void
NonLocalMeansImageFilter <TInputImage>
::BeforeThreadedGenerateData()
{
    Superclass::BeforeThreadedGenerateData();

    this->computeAverageLocalVariance();
    this->computeMeanAndVarImages();
    m_maxAbsDisp = std::floor((double)(m_SearchNeighborhood / m_SearchStepSize)) * m_SearchStepSize;
}

template <class TInputImage>
void
NonLocalMeansImageFilter <TInputImage>
::computeMeanAndVarImages()
{
    typedef MeanAndVarianceImagesFilter<InputImageType, OutputImageType>
            MeanAndVarianceImagesFilterType;
    typename MeanAndVarianceImagesFilterType::Pointer filter = MeanAndVarianceImagesFilterType::New();
    filter->SetInput(this->GetInput());
    typename InputImageType::SizeType radius;

    for (unsigned int j = 0;j < InputImageDimension;++j)
    {
        radius[j] = m_PatchHalfSize;
    }

    filter->SetRadius(radius);
    filter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    filter->Update();
    m_meanImage = filter->GetMeanImage();
    m_varImage = filter->GetVarImage();
}

template <class TInputImage>
void
NonLocalMeansImageFilter <TInputImage>
::computeAverageLocalVariance()
{
    typedef itk::ImageRegionConstIteratorWithIndex< InputImageType > InIteratorType;
    InputImageRegionType largestRegion (this->GetInput()->GetLargestPossibleRegion());
    InIteratorType  dataIterator (this->GetInput(), largestRegion);

    double averageLocalSignal, diffSignal;
    double baseSignal;
    typename InputImageRegionType::IndexType baseIndex;

    double averageCovariance = 0;
    unsigned int numEstimations = 0;
    unsigned int numLocalPixels = 2 * InputImageDimension;

    while (!dataIterator.IsAtEnd())
    {
        baseSignal = static_cast<double>(dataIterator.Get());
        baseIndex = dataIterator.GetIndex();
        averageLocalSignal = 0;

        typename InputImageRegionType::IndexType valueIndex;
        for (unsigned int d = 0; d < InputImageDimension; ++d)
        {
            valueIndex = baseIndex;
            int tmpIndex = baseIndex[d] - m_localNeighborhood;
            valueIndex[d] = std::max(tmpIndex,0);
            averageLocalSignal += static_cast<double> (this->GetInput()->GetPixel(valueIndex));

            valueIndex = baseIndex;
            tmpIndex = baseIndex[d] + m_localNeighborhood;
            int maxIndex = largestRegion.GetSize()[d] - 1;
            valueIndex[d] = std::min(tmpIndex, maxIndex);
            averageLocalSignal += static_cast<double> (this->GetInput()->GetPixel(valueIndex));
        }

        averageLocalSignal /= numLocalPixels;
        diffSignal = sqrt(numLocalPixels / (numLocalPixels + 1.0)) * (baseSignal - averageLocalSignal);

        averageCovariance += diffSignal * diffSignal;

        ++numEstimations;
        ++dataIterator;
    }

    // Now divide by number of estimations and compute average variance
    m_noiseCovariance = averageCovariance / numEstimations;
}

template < class TInputImage>
void
NonLocalMeansImageFilter < TInputImage >
::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    // Allocate output
    typename OutputImageType::Pointer output = this->GetOutput();
    typename InputImageType::Pointer input = const_cast<InputImageType *> (this->GetInput());

    typedef itk::ImageRegionConstIterator< InputImageType > InIteratorType;
    typedef itk::ImageRegionIteratorWithIndex< OutputImageType > OutRegionIteratorType;

    InIteratorType inputIterator(input, outputRegionForThread);
    OutRegionIteratorType outputIterator(output, outputRegionForThread);

    std::vector <InputPixelType> databaseSamples;
    std::vector <double> databaseWeights;

    typedef anima::NonLocalMeansPatchSearcher <TInputImage, OutputImageType> PatchSearcherType;

    PatchSearcherType patchSearcher;
    patchSearcher.SetPatchHalfSize(m_PatchHalfSize);
    patchSearcher.SetSearchStepSize(m_SearchStepSize);
    patchSearcher.SetMaxAbsDisp(m_maxAbsDisp);
    patchSearcher.SetInputImage(input);
    patchSearcher.SetBetaParameter(m_BetaParameter);
    patchSearcher.SetNoiseCovariance(m_noiseCovariance);
    patchSearcher.SetWeightThreshold(m_WeightThreshold);
    patchSearcher.SetMeanImage(m_meanImage);
    patchSearcher.SetVarImage(m_varImage);
    patchSearcher.SetMeanMinThreshold(m_MeanMinThreshold);
    patchSearcher.SetVarMinThreshold(m_VarMinThreshold);

    while (!outputIterator.IsAtEnd())
    {
        patchSearcher.UpdateAtPosition(outputIterator.GetIndex());

        databaseSamples = patchSearcher.GetDatabaseSamples();
        databaseWeights = patchSearcher.GetDatabaseWeights();

        //Compute weighted mean of databaseSamples
        double average = 0, sum = 0, w_max = 0;

        switch (m_WeightMethod)
        {
        case EXP:
            for (unsigned int d = 0;d < databaseSamples.size();++d)
            {
                average += databaseSamples[d] * databaseWeights[d];
                sum += databaseWeights[d];

                if (w_max < databaseWeights[d])
                    w_max = databaseWeights[d];
            }

            if (sum != 0)
                outputIterator.Set((average + w_max * inputIterator.Get()) / (sum + w_max));
            else
                outputIterator.Set(inputIterator.Get());

            break;

        case RICIAN:
            for (unsigned int d=0; d < databaseSamples.size(); d++)
            {
                average += databaseWeights[d] * (databaseSamples[d] * databaseSamples[d]);
                sum += databaseWeights[d];
                if (w_max < databaseWeights[d])
                    w_max = databaseWeights[d];
            }

            if (sum != 0)
            {
                double t = ((average + (inputIterator.Get() * inputIterator.Get())
                             * w_max) / (sum + w_max)) - (2.0 * m_noiseCovariance);

                if (t < 0)
                    t = 0;

                outputIterator.Set(std::sqrt(t));
            }
            else
                outputIterator.Set(inputIterator.Get());

            break;
        }

        this->IncrementNumberOfProcessedPoints();
        ++outputIterator;
        ++inputIterator;
    }
}

} // end of namespace anima
