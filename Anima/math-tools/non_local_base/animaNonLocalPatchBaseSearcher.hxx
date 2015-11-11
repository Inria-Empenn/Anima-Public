#pragma once
#include "animaNonLocalPatchBaseSearcher.h"

#include <itkImageRegionConstIteratorWithIndex.h>

namespace anima
{

template <class ImageType>
NonLocalPatchBaseSearcher <ImageType>
::NonLocalPatchBaseSearcher()
{
    m_PatchHalfSize = 1;
    m_SearchStepSize = 1;
    m_MaxAbsDisp = 3;
    m_WeightThreshold = 0.0;

    m_InputImage = 0;
}

template <class ImageType>
void
NonLocalPatchBaseSearcher <ImageType>
::UpdateAtPosition(const IndexType &dataIndex)
{
    m_DatabaseWeights.clear();
    m_DatabaseSamples.clear();

    IndexType blockIndex, dispIndex, movingIndex;
    SizeType blockSize, dispSize;
    ImageRegionType largestImageRegion = this->GetInputImage()->GetLargestPossibleRegion();
    ImageRegionType blockRegionMoving;

    for (unsigned int d = 0; d < ImageType::ImageDimension;++d)
    {
        int tmpIndex = dataIndex[d] - m_PatchHalfSize;
        blockIndex[d] = std::max(0, tmpIndex);

        int maxSize = largestImageRegion.GetSize()[d] - 1;
        int tmpSize = dataIndex[d] + m_PatchHalfSize;
        blockSize[d] = std::min(maxSize, tmpSize) - blockIndex[d] + 1;

        tmpIndex = dataIndex[d] - m_MaxAbsDisp;
        dispIndex[d] = std::max(0, tmpIndex);

        tmpSize = dataIndex[d] + m_MaxAbsDisp;
        dispSize[d] = std::min(maxSize,tmpSize) - dispIndex[d] + 1;
    }

    ImageRegionType blockRegion;
    blockRegion.SetIndex(blockIndex);
    blockRegion.SetSize(blockSize);

    ImageRegionType dispRegion;
    dispRegion.SetIndex(dispIndex);
    dispRegion.SetSize(dispSize);

    typedef itk::ImageRegionConstIteratorWithIndex <ImageType> InIteratorType;
    InIteratorType dispIt(m_InputImage, dispRegion);

    IndexType dispBaseIndex = dispIt.GetIndex();
    IndexType dispCurIndex;
    while (!dispIt.IsAtEnd())
    {
        dispCurIndex = dispIt.GetIndex();

        bool onSearchStepSize(true), movingRegionIsValid(true), isCentralIndex(true);
        for (unsigned int d = 0; d < ImageType::ImageDimension && onSearchStepSize && movingRegionIsValid; ++d)
        {
            //if the iterator isn't at a search step we won't do anything
            if ((dispCurIndex[d] - dispBaseIndex[d]) % m_SearchStepSize)
            {
                onSearchStepSize = false;
                break;
            }
            else
            {
                movingIndex[d] =  blockIndex[d] + (dispCurIndex[d] - dataIndex[d]);
                unsigned int maxBlock = movingIndex[d] + blockSize[d];
                //if movingRegion overfill largestRegion, we won't compute it
                if (maxBlock > largestImageRegion.GetSize()[d] || movingIndex[d] < 0)
                {
                    movingRegionIsValid = false;
                    break;
                }
            }

            if (dispCurIndex[d] != dataIndex[d])
                isCentralIndex = false;
        }

        if (movingRegionIsValid && onSearchStepSize && (!isCentralIndex))
        {
            this->UpdateAdditionalProperties(blockRegionMoving);
            // Should we compute the weight value of this patch ?
            if (this->testPatchConformity(dataIndex,dispCurIndex))
            {
                blockRegionMoving.SetIndex(movingIndex);
                blockRegionMoving.SetSize(blockRegion.GetSize());

                double weightValue = this->computeWeightValue(blockRegion, blockRegionMoving);
                if (weightValue > m_WeightThreshold)
                {
                    m_DatabaseWeights.push_back(weightValue);
                    // Getting center index value
                    m_DatabaseSamples.push_back(dispIt.Get());
                }
            }
        }

        ++dispIt;
    }
}

} // end namespace anima
