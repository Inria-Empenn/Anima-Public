#pragma once
#include "animaLocalPatchMeanDistanceImageFilter.h"

#include <itkImageRegionIteratorWithIndex.h>
#include <itkSymmetricEigenAnalysis.h>

#include <animaVectorImagePatchStatistics.h>

namespace anima
{

template <class PixelScalarType>
void
LocalPatchMeanDistanceImageFilter<PixelScalarType>
::BeforeThreadedGenerateData ()
{
    Superclass::BeforeThreadedGenerateData();

    // Checking consistency of the data and parameters

    unsigned int nbInputs = this->GetNumberOfIndexedInputs();
    if (nbInputs <= 1)
        itkExceptionMacro("Error: Not enough inputs available... Exiting..." );
}

template <class PixelScalarType>
void
LocalPatchMeanDistanceImageFilter<PixelScalarType>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionIteratorWithIndex< OutputImageType > OutRegionIteratorType;
    typedef itk::ImageRegionConstIteratorWithIndex < MaskImageType > MaskRegionIteratorType;

    OutRegionIteratorType outMeanIterator(this->GetOutput(0), outputRegionForThread);
    OutRegionIteratorType outStdIterator(this->GetOutput(1), outputRegionForThread);

    MaskRegionIteratorType maskIterator (this->GetComputationMask(), outputRegionForThread);

    unsigned int numSamplesDatabase = this->GetNumberOfIndexedInputs();

    unsigned int ndim = this->GetInput(0)->GetNumberOfComponentsPerPixel();
    CovarianceType baseVar(ndim,ndim,0);
    VectorType patchMean(ndim);

    std::vector <unsigned int> numPixels(numSamplesDatabase, 0);
    std::vector <VectorType> meanVectors(numSamplesDatabase, patchMean);
    std::vector <CovarianceType> varianceVector(numSamplesDatabase, baseVar);
    OutputImageRegionType tmpBlockRegion;

    unsigned int numDistances = numSamplesDatabase * (numSamplesDatabase + 1) / 2 - numSamplesDatabase;

    InputImageIndexType curIndex;
    OutputImageRegionType largestRegionOut = this->GetOutput(0)->GetLargestPossibleRegion();

    while (!maskIterator.IsAtEnd())
    {
        if (maskIterator.Get() == 0)
        {
            outMeanIterator.Set(0.0);
            outStdIterator.Set(0.0);

            ++outMeanIterator;
            ++outStdIterator;
            ++maskIterator;
            continue;
        }

        curIndex = maskIterator.GetIndex();

        for (unsigned int i = 0;i < 3;++i)
        {
            tmpBlockRegion.SetIndex(i,std::max(0,(int)curIndex[i] - (int)m_PatchHalfSize));
            tmpBlockRegion.SetSize(i,std::min((unsigned int)(largestRegionOut.GetSize()[i] - 1),(unsigned int)(curIndex[i] + m_PatchHalfSize)) - tmpBlockRegion.GetIndex(i) + 1);
        }

        for (unsigned int i = 0;i < numSamplesDatabase;++i)
            numPixels[i] = anima::computePatchMeanAndCovariance(this->GetInput(i),tmpBlockRegion,meanVectors[i],varianceVector[i]);

        double meanDist = 0;
        double varDist = 0;
        for (unsigned int i = 0;i < numSamplesDatabase;++i)
            for (unsigned int j = i+1;j < numSamplesDatabase;++j)
            {
                double tmpDist = anima::VectorMeansTest(meanVectors[i], meanVectors[j], numPixels[i], numPixels[j],
                                                        varianceVector[i], varianceVector[j]);
                meanDist += tmpDist;
                varDist += tmpDist * tmpDist;
            }

        varDist /= numDistances;
        meanDist /= numDistances;
        varDist -= meanDist * meanDist;
        varDist *= numDistances / (numDistances - 1.0);

        outMeanIterator.Set(meanDist);
        outStdIterator.Set(std::sqrt(varDist));

        ++outMeanIterator;
        ++outStdIterator;
        ++maskIterator;
    }
}

} // end namespace anima
