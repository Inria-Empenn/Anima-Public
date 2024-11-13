#pragma once
#include "animaBlockMatchInitializer.h"

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIteratorWithIndex.h>

#include <itkExpNegativeImageFilter.h>
#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkPoolMultiThreader.h>

namespace anima
{

template <class PixelType, unsigned int NDimensions>
std::vector < typename BlockMatchingInitializer<PixelType,NDimensions>::ImageRegionType > &
BlockMatchingInitializer<PixelType,NDimensions>
::GetOutput()
{
    this->Update();

    return m_Output;
}

template <class PixelType, unsigned int NDimensions>
std::vector < typename BlockMatchingInitializer<PixelType,NDimensions>::PointType > &
BlockMatchingInitializer<PixelType,NDimensions>
::GetOutputPositions()
{
    this->Update();

    return m_OutputPositions;
}

template <class PixelType, unsigned int NDimensions>
std::vector <unsigned int> &
BlockMatchingInitializer<PixelType,NDimensions>::GetMaskStartingIndexes()
{
    if (m_MaskStartingIndexes.size() == 0)
        this->Update();

    return m_MaskStartingIndexes;
}

template <class PixelType, unsigned int NDimensions>
void
BlockMatchingInitializer<PixelType,NDimensions>
::AddReferenceImage(itk::ImageBase <NDimensions> *refImage)
{
    ScalarImageType * tmpPtr = dynamic_cast <ScalarImageType *> (refImage);

    if (tmpPtr)
    {
        m_ReferenceScalarImages.push_back(tmpPtr);
        m_ReferenceScalarImages.back()->DisconnectPipeline();
    }
    else
    {
        m_ReferenceVectorImages.push_back(dynamic_cast <VectorImageType *> (refImage));
        m_ReferenceVectorImages.back()->DisconnectPipeline();
    }
}

template <class PixelType, unsigned int NDimensions>
void
BlockMatchingInitializer<PixelType,NDimensions>
::AddGenerationMask(MaskImageType *mask)
{
    if (mask)
        m_GenerationMasks.push_back(mask);
}

template <class PixelType, unsigned int NDimensions>
itk::ImageBase <NDimensions> *
BlockMatchingInitializer<PixelType,NDimensions>
::GetFirstReferenceImage()
{
    if (m_ReferenceScalarImages.size() != 0)
        return m_ReferenceScalarImages[0];

    if (m_ReferenceVectorImages.size() != 0)
        return m_ReferenceVectorImages[0];

    return NULL;
}

template <class PixelType, unsigned int NDimensions>
void BlockMatchingInitializer<PixelType,NDimensions>
::SetPercentageKept(double val)
{
    if (val != m_PercentageKept)
    {
        m_PercentageKept = val;
        m_UpToDate = false;
    }
}

template <class PixelType, unsigned int NDimensions>
void BlockMatchingInitializer<PixelType,NDimensions>
::SetScalarVarianceThreshold(double val)
{
    if (val != m_ScalarVarianceThreshold)
    {
        m_ScalarVarianceThreshold = val;
        m_UpToDate = false;
    }
}

template <class PixelType, unsigned int NDimensions>
void BlockMatchingInitializer<PixelType,NDimensions>
::SetOrientedModelVarianceThreshold(double val)
{
    if (val != m_OrientedModelVarianceThreshold)
    {
        m_OrientedModelVarianceThreshold = val;
        m_UpToDate = false;
    }
}

template <class PixelType, unsigned int NDimensions>
void BlockMatchingInitializer<PixelType,NDimensions>
::SetRequestedRegion(const ImageRegionType &val)
{
    if (val != m_RequestedRegion)
    {
        m_RequestedRegion = val;
        m_UpToDate = false;
    }
}

template <class PixelType, unsigned int NDimensions>
void BlockMatchingInitializer<PixelType,NDimensions>
::SetBlockSpacing(unsigned int val)
{
    if (val != m_BlockSpacing)
    {
        m_BlockSpacing = val;
        m_UpToDate = false;
    }
}

template <class PixelType, unsigned int NDimensions>
void BlockMatchingInitializer<PixelType,NDimensions>
::SetBlockSize(unsigned int val)
{
    if (val != m_BlockSize)
    {
        m_BlockSize = val;
        m_UpToDate = false;
    }
}

template <class PixelType, unsigned int NDimensions>
void BlockMatchingInitializer<PixelType,NDimensions>
::Update()
{
    if (m_GenerationMasks.size() == 0)
    {
        // Not taking origin, orientation, we don't care about it
        MaskImagePointer wholeMask = MaskImageType::New();
        wholeMask->Initialize();
        wholeMask->SetRegions(m_RequestedRegion);
        wholeMask->Allocate();
        wholeMask->FillBuffer(1);

        m_GenerationMasks.push_back(wholeMask);
        m_UpToDate = false;
    }

    if (m_UpToDate)
        return;

    m_Output.clear();
    m_OutputPositions.clear();
    m_MaskStartingIndexes.resize(m_GenerationMasks.size());
    for (unsigned int i = 0;i < m_GenerationMasks.size();++i)
    {
        m_MaskStartingIndexes[i] = m_Output.size();
        this->ComputeBlocksOnGenerationMask(i);
    }

    m_UpToDate = true;
}

template <class PixelType, unsigned int NDimensions>
void
BlockMatchingInitializer<PixelType,NDimensions>
::ComputeBlocksOnGenerationMask(unsigned int maskIndex)
{
    itk::PoolMultiThreader::Pointer threaderBlockGenerator = itk::PoolMultiThreader::New();

    BlockGeneratorThreadStruct *tmpStr = 0;
    this->InitializeThreading(maskIndex,tmpStr);

    threaderBlockGenerator->SetNumberOfWorkUnits(this->GetNumberOfThreads());
    threaderBlockGenerator->SetSingleMethod(this->ThreadBlockGenerator,tmpStr);
    threaderBlockGenerator->SingleMethodExecute();

    m_Output.clear();
    m_OutputPositions.clear();
    unsigned int totalNumberOfBlocks = 0;
    unsigned int realNumberOfBlocks = 0;

    for (unsigned int i = 0;i < this->GetNumberOfThreads();++i)
    {
        realNumberOfBlocks += tmpStr->tmpOutput[i].size();
        totalNumberOfBlocks += tmpStr->totalNumberOfBlocks[i];
    }

    double percentageBlocksKept = (double) realNumberOfBlocks / totalNumberOfBlocks;

    if (percentageBlocksKept > m_PercentageKept)
    {
        std::vector < std::pair <double, std::pair <PointType, ImageRegionType> > > sortVector;
        for (unsigned int i = 0;i < this->GetNumberOfThreads();++i)
            for (unsigned int j = 0;j < tmpStr->tmpOutput[i].size();++j)
            {
                std::pair <PointType, ImageRegionType> tmpPair(tmpStr->blocks_positions[i][j],tmpStr->tmpOutput[i][j]);
                sortVector.push_back(std::make_pair(tmpStr->blocks_variances[i][j],tmpPair));
            }

        unsigned int numRemoved = std::floor((1.0 - m_PercentageKept) * totalNumberOfBlocks);
        std::partial_sort(sortVector.begin(),sortVector.begin() + numRemoved,sortVector.end(),pair_comparator());

        for (unsigned int i = numRemoved;i < sortVector.size();++i)
        {
            m_Output.push_back(sortVector[i].second.second);
            m_OutputPositions.push_back(sortVector[i].second.first);
        }
    }
    else
    {
        for (unsigned int i = 0;i < this->GetNumberOfThreads();++i)
        {
            m_Output.insert(m_Output.end(), tmpStr->tmpOutput[i].begin(), tmpStr->tmpOutput[i].end());
            m_OutputPositions.insert(m_OutputPositions.end(), tmpStr->blocks_positions[i].begin(),
                                     tmpStr->blocks_positions[i].end());
        }
    }

    delete tmpStr;
}

template <class PixelType, unsigned int NDimensions>
void
BlockMatchingInitializer<PixelType, NDimensions>
::InitializeThreading(unsigned int maskIndex, BlockGeneratorThreadStruct *&workStr)
{
    ImageRegionType workRegion;
    IndexType minIndex, maxIndex, tmpIndex;
    minIndex = m_RequestedRegion.GetIndex() + m_RequestedRegion.GetSize();
    maxIndex = m_RequestedRegion.GetIndex();
    typedef itk::ImageRegionIteratorWithIndex <MaskImageType> MaskIteratorType;
    MaskIteratorType maskItr(m_GenerationMasks[maskIndex],m_RequestedRegion);

    while (!maskItr.IsAtEnd())
    {
        if (maskItr.Get() != 0)
        {
            tmpIndex = maskItr.GetIndex();
            for (unsigned int i = 0;i < NDimensions;++i)
            {
                if (tmpIndex[i] < minIndex[i])
                    minIndex[i] = tmpIndex[i];

                if (tmpIndex[i] > maxIndex[i])
                    maxIndex[i] = tmpIndex[i];
            }
        }

        ++maskItr;
    }

    for (unsigned int i = 0;i < NDimensions;++i)
    {
        workRegion.SetIndex(i,minIndex[i]);
        workRegion.SetSize(i,maxIndex[i] - minIndex[i] + 1);
    }

    if (workStr == 0)
        workStr = new BlockGeneratorThreadStruct;

    std::vector <unsigned int> totalNbBlocks(NDimensions);
    workStr->blockStartOffsets.resize(NDimensions);

    for (unsigned int i = 0;i < NDimensions;++i)
    {
        totalNbBlocks[i] = std::floor((double)(workRegion.GetSize()[i] / this->GetBlockSpacing()) + 0.5);
        if (totalNbBlocks[i] < 1)
            totalNbBlocks[i] = 1;

        unsigned int spaceRequired = workRegion.GetSize()[i] - (totalNbBlocks[i] - 1) * this->GetBlockSpacing() - 1;

        workStr->blockStartOffsets[i] = workRegion.GetIndex()[i] + std::floor(spaceRequired / 2.0);
    }

    unsigned int nb_blocks_per_thread = (unsigned int) std::floor((double)(totalNbBlocks[NDimensions-1] / this->GetNumberOfThreads()));
    if (nb_blocks_per_thread < 1)
    {
        nb_blocks_per_thread = 1;
        this->SetNumberOfThreads(totalNbBlocks[NDimensions-1]);
    }

    workStr->startBlocks.resize(this->GetNumberOfThreads());
    workStr->nb_blocks.resize(this->GetNumberOfThreads());
    workStr->maskIndex = maskIndex;

    unsigned int currentCount = 0;
    for (unsigned int i = 0;i < this->GetNumberOfThreads();++i)
    {
        std::vector <unsigned int> startBlock(NDimensions,0);
        startBlock[NDimensions-1] = currentCount;
        workStr->startBlocks[i] = startBlock;

        std::vector <unsigned int> numBlockVec(NDimensions);
        for (unsigned int j = 0;j < NDimensions-1;++j)
            numBlockVec[j] = totalNbBlocks[j];

        unsigned int numBlocksForThread = nb_blocks_per_thread;
        if (numBlocksForThread + currentCount > totalNbBlocks[NDimensions-1])
            numBlocksForThread = totalNbBlocks[NDimensions-1] - currentCount;

        //We may be off by a few planes of blocks at the end
        if (i == this->GetNumberOfThreads()-1)
            numBlocksForThread = totalNbBlocks[NDimensions-1] - currentCount;

        numBlockVec[NDimensions-1] = numBlocksForThread;
        workStr->nb_blocks[i] = numBlockVec;

        currentCount += numBlocksForThread;
    }

    workStr->Filter = this;
    workStr->tmpOutput.resize(this->GetNumberOfThreads());
    workStr->totalNumberOfBlocks.resize(this->GetNumberOfThreads());
    workStr->blocks_positions.resize(this->GetNumberOfThreads());
    workStr->blocks_variances.resize(this->GetNumberOfThreads());

    for (unsigned int i = 0;i < this->GetNumberOfThreads();++i)
    {
        workStr->tmpOutput[i].clear();
        workStr->totalNumberOfBlocks[i] = 0;
        workStr->blocks_positions[i].clear();
        workStr->blocks_variances[i].clear();
    }
}

template <class PixelType, unsigned int NDimensions>
ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION
BlockMatchingInitializer<PixelType,NDimensions>
::ThreadBlockGenerator(void *arg)
{
    itk::MultiThreaderBase::WorkUnitInfo *threadArgs = (itk::MultiThreaderBase::WorkUnitInfo *)arg;

    unsigned int nbThread = threadArgs->WorkUnitID;

    BlockGeneratorThreadStruct *tmpStr = (BlockGeneratorThreadStruct *)threadArgs->UserData;

    tmpStr->Filter->RegionBlockGenerator(tmpStr,nbThread);

    return ITK_THREAD_RETURN_DEFAULT_VALUE;
}

template <class PixelType, unsigned int NDimensions>
void BlockMatchingInitializer<PixelType,NDimensions>
::RegionBlockGenerator(BlockGeneratorThreadStruct *workStr, unsigned int threadId)
{
    typename ImageRegionType::IndexType startPosition, blockPosition;
    typename ImageRegionType::SizeType blockSize;
    PointType blockOrigin;
    ImageRegionType largestRegion = this->GetFirstReferenceImage()->GetLargestPossibleRegion();

    ImageRegionType tmpBlock;
    int indexPos;
    double blockVar = 0;
    workStr->totalNumberOfBlocks[threadId] = 0;
    workStr->blocks_variances[threadId].clear();
    workStr->blocks_positions[threadId].clear();

    unsigned int block_half_size = std::floor ((this->GetBlockSize() - 1) / 2.0);

    bool continueLoop = true;
    std::vector <unsigned int> positionCounter(NDimensions,0);

    while (continueLoop)
    {
        for (unsigned int i = 0;i < NDimensions;++i)
        {
            blockSize[i] = this->GetBlockSize();
            indexPos = workStr->blockStartOffsets[i] + (workStr->startBlocks[threadId][i] + positionCounter[i])*this->GetBlockSpacing() - block_half_size;
            if (indexPos < largestRegion.GetIndex()[i])
            {
                blockSize[i] += indexPos - largestRegion.GetIndex()[i];
                indexPos = largestRegion.GetIndex()[i];
            }
            startPosition[i] = indexPos;
            if (startPosition[i] + blockSize[i] > largestRegion.GetIndex()[i] + largestRegion.GetSize()[i])
                blockSize[i] = largestRegion.GetIndex()[i] + largestRegion.GetSize()[i] - startPosition[i];
        }

        tmpBlock.SetIndex(startPosition);
        tmpBlock.SetSize(blockSize);
        workStr->totalNumberOfBlocks[threadId]++;

        if (this->CheckBlockConditions(tmpBlock,blockVar,workStr,threadId))
        {
            for (unsigned int i = 0;i < NDimensions;++i)
                blockPosition[i] = workStr->blockStartOffsets[i] + (workStr->startBlocks[threadId][i] + positionCounter[i])*this->GetBlockSpacing();

            if (m_GenerationMasks[workStr->maskIndex]->GetPixel(blockPosition) == 0)
            {
                continueLoop = this->ProgressCounter(positionCounter,workStr->nb_blocks[threadId]);
                continue;
            }

            workStr->tmpOutput[threadId].push_back(tmpBlock);
            workStr->blocks_variances[threadId].push_back(blockVar);

            this->GetFirstReferenceImage()->TransformIndexToPhysicalPoint(blockPosition,blockOrigin);
            workStr->blocks_positions[threadId].push_back(blockOrigin);
        }

        continueLoop = this->ProgressCounter(positionCounter,workStr->nb_blocks[threadId]);
    }
}

template <class PixelType, unsigned int NDimensions>
bool
BlockMatchingInitializer<PixelType,NDimensions>
::ProgressCounter(std::vector <unsigned int> &counter, std::vector <unsigned int> & bounds)
{
    unsigned int pos = counter.size() - 1;

    bool needUpdate = false;
    do
    {
        counter[pos]++;
        needUpdate = false;

        if ((counter[pos] == bounds[pos])&&(pos > 0))
        {
            counter[pos] = 0;
            --pos;
            needUpdate = true;
        }
    } while (needUpdate);

    if (counter[pos] >= bounds[pos])
        return false;

    return true;
}

template <class PixelType, unsigned int NDimensions>
bool
BlockMatchingInitializer<PixelType,NDimensions>
::CheckBlockConditions(ImageRegionType &region, double &blockVariance, BlockGeneratorThreadStruct *workStr,
                       unsigned int threadId)
{
    ImageRegionType refRegion = this->GetFirstReferenceImage()->GetLargestPossibleRegion();

    for (unsigned int i = 0;i < NDimensions;++i)
    {
        if (refRegion.GetIndex()[i] > region.GetIndex()[i])
            return false;

        if (refRegion.GetIndex()[i] + refRegion.GetSize()[i] < region.GetIndex()[i] + region.GetSize()[i])
            return false;
    }

    blockVariance = 0;
    double tmpVar = 0;

    for (unsigned int i = 0;i < m_ReferenceScalarImages.size();++i)
    {
        if (!this->CheckScalarVariance(m_ReferenceScalarImages[i],region,tmpVar))
            return false;

        if (tmpVar > blockVariance)
            blockVariance = tmpVar;
    }

    for (unsigned int i = 0;i < m_ReferenceVectorImages.size();++i)
    {
        if (!this->CheckOrientedModelVariance(i,region,tmpVar,workStr,threadId))
            return false;

        if (tmpVar > blockVariance)
            blockVariance = tmpVar;
    }

    // Reaches here if block respects criterions on all images
    return true;
}

template <class PixelType, unsigned int NDimensions>
bool BlockMatchingInitializer<PixelType,NDimensions>
::CheckOrientedModelVariance(unsigned int imageIndex, ImageRegionType &region, double &blockVariance,
                             BlockGeneratorThreadStruct *workStr, unsigned int threadId)
{
    VectorImageType *refImage = m_ReferenceVectorImages[imageIndex];
    itk::ImageRegionConstIterator <VectorImageType> refItr(refImage,region);
    typedef typename VectorImageType::PixelType VectorType;

    unsigned int nbPts = 0;
    unsigned int vectorSize = refImage->GetNumberOfComponentsPerPixel();

    VectorType meanVal(vectorSize);
    meanVal.Fill(0);
    std::vector <double> blockCovariance(vectorSize, 0.0);

    while (!refItr.IsAtEnd())
    {
        VectorType tmpVal = refItr.Get();
        meanVal += tmpVal;

        ++nbPts;
        ++refItr;
    }

    meanVal /= nbPts;

    refItr.GoToBegin();

    while (!refItr.IsAtEnd())
    {
        VectorType tmpVal = refItr.Get();

        for (unsigned int j = 0;j < vectorSize;++j)
            blockCovariance[j] += (tmpVal[j] - meanVal[j]) * (tmpVal[j] - meanVal[j]);

        ++refItr;
    }

    blockVariance = 0;

    if (nbPts <= 1)
        return false;

    for (unsigned int j = 0;j < vectorSize;++j)
    {
        double tmp = blockCovariance[j] / (nbPts - 1.0);
        blockVariance += tmp;
    }

    blockVariance /= vectorSize;

    if (blockVariance > this->GetOrientedModelVarianceThreshold())
        return true;
    else
        return false;
}

template <class PixelType, unsigned int NDimensions>
bool BlockMatchingInitializer<PixelType,NDimensions>
::CheckScalarVariance(ScalarImageType *refImage, ImageRegionType &region, double &blockVariance)
{
    itk::ImageRegionConstIterator <ScalarImageType> refItr(refImage,region);
    blockVariance = 0;
    double meanVal = 0;

    unsigned int nbPts = 0;
    while (!refItr.IsAtEnd())
    {
        double tmpVal = refItr.Get();
        meanVal += tmpVal;

        ++nbPts;
        ++refItr;
    }

    if (nbPts <= 1)
        return false;

    meanVal /= nbPts;
    refItr.GoToBegin();
    while (!refItr.IsAtEnd())
    {
        double tmpVal = refItr.Get();
        blockVariance += (meanVal - tmpVal) * (meanVal - tmpVal);

        ++refItr;
    }

    blockVariance /= (nbPts - 1.0);

    if (blockVariance > this->GetScalarVarianceThreshold())
        return true;
    else
        return false;
}

}// end of namespace anima
