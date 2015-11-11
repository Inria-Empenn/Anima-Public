#pragma once
#include "animaNonLocalMeansPatchSearcher.h"

namespace anima
{

template <class ImageType, class DataImageType>
NonLocalMeansPatchSearcher <ImageType, DataImageType>
::NonLocalMeansPatchSearcher()
{
    m_BetaParameter = 1.0;
    m_NoiseCovariance = 1.0;

    m_MeanMinThreshold = 0.95;
    m_VarMinThreshold = 0.5;
}

template <class ImageType, class DataImageType>
bool
NonLocalMeansPatchSearcher <ImageType, DataImageType>
::TestPatchConformity(unsigned int index, const IndexType &refIndex, const IndexType &movingIndex)
{
    double refMeanValue = m_MeanImage->GetPixel(refIndex);
    double floMeanValue = m_MeanImage->GetPixel(movingIndex);

    double refVarValue = m_VarImage->GetPixel(refIndex);
    double floVarValue = m_VarImage->GetPixel(movingIndex);

    double meanRate = refMeanValue / floMeanValue;
    double varianceRate = refVarValue / floVarValue;

    // Should we compute the weight value of this patch ?
    if ( ( meanRate > m_MeanMinThreshold ) && ( meanRate < ( 1.0 / m_MeanMinThreshold ) ) &&
         ( varianceRate > m_VarMinThreshold ) && ( varianceRate < ( 1.0 / m_VarMinThreshold ) ) )
        return true;

    return false;
}

template <class ImageType, class DataImageType>
double
NonLocalMeansPatchSearcher <ImageType, DataImageType>
::ComputeWeightValue(unsigned int index, ImageRegionType &refPatch, ImageRegionType &movingPatch)
{
    typedef itk::ImageRegionConstIteratorWithIndex< ImageType > InIteratorType;

    InIteratorType tmpIt (this->GetInputImage(), refPatch);
    InIteratorType tmpMovingIt (this->GetComparisonImage(index), movingPatch);

    double tmpDiffValue;

    double weightValue = 0.0;
    unsigned int numVoxels = 0;

    while (!tmpIt.IsAtEnd())
    {
        tmpDiffValue = (double)tmpIt.Get() - (double)tmpMovingIt.Get();
        weightValue += tmpDiffValue * tmpDiffValue;

        ++numVoxels;
        ++tmpIt;
        ++tmpMovingIt;
    }

    weightValue = std::exp(- weightValue / (2.0 * m_BetaParameter * m_NoiseCovariance * numVoxels));
    return weightValue;
}

} // end namespace anima
