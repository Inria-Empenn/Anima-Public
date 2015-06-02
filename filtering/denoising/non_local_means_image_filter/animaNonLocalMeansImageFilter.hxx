#pragma once
#include "animaNonLocalMeansImageFilter.h"

#include <animaMeanAndVarianceImagesFilter.h>

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
    m_maxAbsDisp = floor((float)(m_SearchNeighborhood / m_SearchStepSize)) * m_SearchStepSize;
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

    for (unsigned int j = 0 ; (j < InputImageDimension); ++j)
    {
        radius[j] = m_PatchHalfSize;
    }

    filter->SetRadius(radius);
    filter->SetNumberOfThreads(this->GetNumberOfThreads());
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
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{

    // support progress methods/callbacks
    itk::ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
     // Allocate output
    typename OutputImageType::Pointer output = this->GetOutput();
    typename  InputImageType::ConstPointer input  = this->GetInput();

    typedef itk::ImageRegionConstIteratorWithIndex< InputImageType > InIteratorType;
    typedef itk::ImageRegionIteratorWithIndex< OutputImageType > OutRegionIteratorType;

    InIteratorType inputIterator(input, outputRegionForThread);
    OutRegionIteratorType outputIterator(output, outputRegionForThread);
    OutRegionIteratorType meanIterator(m_meanImage, outputRegionForThread);
    OutRegionIteratorType varIterator(m_varImage, outputRegionForThread);


    std::vector <InputPixelType> databaseSamples;
    std::vector <double> databaseWeights;

    OutputImageRegionType largestRegionOut = output->GetLargestPossibleRegion();
    OutputImageRegionType blockRegion, blockRegionMoving;


    InputImageRegionType dispRegion;
    typename InputImageRegionType::IndexType blockIndex, dispIndex, curIndex, dispCurIndex, movingIndex, dispBaseIndex;
    typename InputImageRegionType::SizeType blockSize, dispSize;

    while (!outputIterator.IsAtEnd())
    {
        databaseWeights.clear();
        databaseSamples.clear();

        curIndex = outputIterator.GetIndex();

        for (unsigned int d = 0; d < InputImageDimension; ++d)
        {
            int tmpIndex = curIndex[d] - m_PatchHalfSize;
            blockIndex[d] = std::max(0, tmpIndex);

            int maxSize = largestRegionOut.GetSize()[d] - 1;
            int tmpSize = curIndex[d] + m_PatchHalfSize;
            blockSize[d] = std::min(maxSize, tmpSize) - blockIndex[d] + 1;

            tmpIndex = curIndex[d] - m_maxAbsDisp;
            dispIndex[d] = std::max(0, tmpIndex);

            tmpSize = curIndex[d] + m_maxAbsDisp;
            dispSize[d] = std::min(maxSize,tmpSize) - dispIndex[d] + 1;
        }

        blockRegion.SetIndex(blockIndex);
        blockRegion.SetSize(blockSize);

        dispRegion.SetIndex(dispIndex);
        dispRegion.SetSize(dispSize);

        InIteratorType  dispIt(input, dispRegion);
        OutRegionIteratorType  dispMeanIt(m_meanImage, dispRegion);
        OutRegionIteratorType  dispVarIt(m_varImage, dispRegion);

        dispBaseIndex = dispIt.GetIndex();
        while (!dispIt.IsAtEnd())
        {
            dispCurIndex = dispIt.GetIndex();

            bool onSearchStepSize(true), movingRegionIsValid(true), isCentralIndex(true);
            for (unsigned int d = 0; d < InputImageDimension && onSearchStepSize && movingRegionIsValid; ++d)
            {
                //if the iterator isn't at a search step we won't do anything
                if ((dispCurIndex[d] - dispBaseIndex[d]) % m_SearchStepSize)
                {
                    onSearchStepSize = false;
                    break;
                }
                else
                {
                    movingIndex[d] =  blockIndex[d] + (dispCurIndex[d] - curIndex[d]);
                    unsigned int maxBlock = movingIndex[d] + blockSize[d];
                    //if movingRegion overfill largestRegion, we won't compute it
                    if (maxBlock > largestRegionOut.GetSize()[d] || movingIndex[d] < 0)
                    {
                        movingRegionIsValid = false;
                        break;
                    }
                }

                if (dispCurIndex[d] != curIndex[d])
                    isCentralIndex = false;
            }

            if (movingRegionIsValid && onSearchStepSize && (!isCentralIndex))
            {
                double meanRate = static_cast<double>(meanIterator.Get()) / static_cast<double>(dispMeanIt.Get());
                double varianceRate = static_cast<double>(varIterator.Get()) / static_cast<double>(dispVarIt.Get());

                // Should we compute the weight value of this patch ?
                if ( ( meanRate > m_MeanMinThreshold ) && ( meanRate < ( 1 / m_MeanMinThreshold ) ) &&
                        ( varianceRate > m_VarMinThreshold ) && ( varianceRate < ( 1 / m_VarMinThreshold ) ) )
                {
                    blockRegionMoving.SetIndex(movingIndex);
                    blockRegionMoving.SetSize(blockRegion.GetSize());

                    double weightValue = this->computeWeightValue(blockRegion, blockRegionMoving);
                    if (weightValue > m_WeightThreshold)
                    {
                        databaseWeights.push_back(weightValue);
                        // Getting center index value
                        databaseSamples.push_back(dispIt.Get());
                    }
                }
            }

            ++dispIt;
            ++dispMeanIt;
            ++dispVarIt;
        }

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

                    outputIterator.Set(sqrt(t));
                }
                else
                    outputIterator.Set(inputIterator.Get());

                break;

        }

        ++outputIterator;
        ++inputIterator;
        ++meanIterator;
        ++varIterator;
        progress.CompletedPixel();
    }

}



template < class TInputImage>
double
NonLocalMeansImageFilter < TInputImage >
::computeWeightValue(const OutputImageRegionType& inputRegion,
                     const OutputImageRegionType& movingRegion)
{
    typedef itk::ImageRegionConstIteratorWithIndex< InputImageType > InIteratorType;

    InIteratorType tmpIt (this->GetInput(), inputRegion);
    InIteratorType tmpMovingIt (this->GetInput(), movingRegion);

    double tmpDiffValue;

    double weightValue = 0.0;
    unsigned int numVoxels = 0;

    tmpIt.GoToBegin();
    tmpMovingIt.GoToBegin();
    while (!tmpIt.IsAtEnd())
    {
        tmpDiffValue = (double)tmpIt.Get() - (double)tmpMovingIt.Get();
        weightValue += tmpDiffValue * tmpDiffValue;

        ++numVoxels;
        ++tmpIt;
        ++tmpMovingIt;
    }
    weightValue = exp(- weightValue / (2.0 * m_BetaParameter * m_noiseCovariance * numVoxels));
    return weightValue;
}

} // end of namespace anima
