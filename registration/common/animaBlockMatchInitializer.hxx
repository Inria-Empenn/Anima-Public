#pragma once

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageFileWriter.h>
#include <itkDanielssonDistanceMapImageFilter.h>

namespace anima
{

template <class PixelType, unsigned int NDimensions>
std::vector < typename BlockMatchingInitializer<PixelType,NDimensions>::ImageRegionType > &
BlockMatchingInitializer<PixelType,NDimensions>::GetOutput()
{
    this->Update();

    return m_Output;
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
void BlockMatchingInitializer<PixelType,NDimensions>::Update()
{
    itk::MultiThreader::Pointer threaderBlockGenerator = itk::MultiThreader::New();

    std::vector <unsigned int> totalNbBlocks(NDimensions);

    for (unsigned int i = 0;i < NDimensions;++i)
        totalNbBlocks[i] = std::max(1,(int) floor((float)(this->GetRequestedRegion().GetSize()[i] / this->GetBlockSpacing())));

    unsigned int nb_blocks_per_thread = (unsigned int) floor((float)(totalNbBlocks[NDimensions-1] / this->GetNumberOfThreads()));
    if (nb_blocks_per_thread < 1)
        nb_blocks_per_thread = 1;

    std::vector < std::vector <unsigned int> > startBlocks(this->GetNumberOfThreads()), nb_blocks(this->GetNumberOfThreads());

    unsigned int currentCount = 0;
    for (unsigned int i = 0;i < this->GetNumberOfThreads();++i)
    {
        std::vector <unsigned int> startBlock(NDimensions,0);
        startBlock[NDimensions-1] = currentCount;
        startBlocks[i] = startBlock;

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
        nb_blocks[i] = numBlockVec;

        currentCount += numBlocksForThread;
    }

    BlockGeneratorThreadStruct *tmpStr = new BlockGeneratorThreadStruct;
    tmpStr->Filter = this;
    tmpStr->startBlocks = startBlocks;
    tmpStr->nb_blocks = nb_blocks;

    tmpStr->tmpOutput.resize(this->GetNumberOfThreads());
    tmpStr->totalNumberOfBlocks.resize(this->GetNumberOfThreads());
    tmpStr->blocks_variances.resize(this->GetNumberOfThreads());

    for (unsigned int i = 0;i < this->GetNumberOfThreads();++i)
    {
        tmpStr->tmpOutput[i].clear();
        tmpStr->totalNumberOfBlocks[i] = 0;
        tmpStr->blocks_variances[i].clear();
    }

    threaderBlockGenerator->SetNumberOfThreads(this->GetNumberOfThreads());
    threaderBlockGenerator->SetSingleMethod(this->ThreadBlockGenerator,tmpStr);
    threaderBlockGenerator->SingleMethodExecute();

    m_Output.clear();
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
        std::vector < std::pair <double, ImageRegionType> > sortVector;
        for (unsigned int i = 0;i < this->GetNumberOfThreads();++i)
            for (unsigned int j = 0;j < tmpStr->tmpOutput[i].size();++j)
                sortVector.push_back(std::make_pair(tmpStr->blocks_variances[i][j],tmpStr->tmpOutput[i][j]));

        unsigned int numRemoved = std::floor((1.0 - m_PercentageKept) * totalNumberOfBlocks);
        std::partial_sort(sortVector.begin(),sortVector.begin() + numRemoved,sortVector.end(),pair_comparator());

        for (unsigned int i = numRemoved;i < sortVector.size();++i)
            m_Output.push_back(sortVector[i].second);
    }
    else
    {
        for (unsigned int i = 0;i < this->GetNumberOfThreads();++i)
            m_Output.insert( m_Output.end(), tmpStr->tmpOutput[i].begin(), tmpStr->tmpOutput[i].end() );
    }

    delete tmpStr;

    if (m_ComputeOuterDam)
        this->ComputeOuterDamFromBlocks();
}

template <class PixelType, unsigned int NDimensions>
ITK_THREAD_RETURN_TYPE
BlockMatchingInitializer<PixelType,NDimensions>
::ThreadBlockGenerator(void *arg)
{
    itk::MultiThreader::ThreadInfoStruct *threadArgs = (itk::MultiThreader::ThreadInfoStruct *)arg;

    unsigned int nbThread = threadArgs->ThreadID;

    BlockGeneratorThreadStruct *tmpStr = (BlockGeneratorThreadStruct *)threadArgs->UserData;

    tmpStr->Filter->RegionBlockGenerator(tmpStr->startBlocks[nbThread],tmpStr->nb_blocks[nbThread],
                                         tmpStr->tmpOutput[nbThread],tmpStr->blocks_variances[nbThread],
                                         tmpStr->totalNumberOfBlocks[nbThread]);

    return NULL;
}

template <class PixelType, unsigned int NDimensions>
void BlockMatchingInitializer<PixelType,NDimensions>::
RegionBlockGenerator(std::vector <unsigned int> &num_start, std::vector <unsigned int> &num_blocks,
                     std::vector <ImageRegionType> &tmpOutput, std::vector <double> &block_variances,
                     unsigned int &total_num_blocks)
{
    typename ImageRegionType::IndexType startPosition;
    typename ImageRegionType::IndexType reqStartIndex = this->GetRequestedRegion().GetIndex();
    typename ImageRegionType::SizeType blockSize;

    for (unsigned int i = 0;i < NDimensions;++i)
    {
        if (this->GetBlockSize() <= this->GetRequestedRegion().GetSize()[i])
            blockSize[i] = this->GetBlockSize();
        else
            blockSize[i] = this->GetRequestedRegion().GetSize()[i];
    }

    ImageRegionType tmpBlock;
    tmpBlock.SetSize(blockSize);
    double blockVar = 0;
    total_num_blocks = 0;
    block_variances.clear();

    if (NDimensions == 2)
    {
        for (unsigned int j = 0;j < num_blocks[1];++j)
        {
            startPosition[1] = reqStartIndex[1] + (num_start[1] + j)*this->GetBlockSpacing();
            for (unsigned int i = 0;i < num_blocks[0];++i)
            {
                startPosition[0] = reqStartIndex[0] + (num_start[0] + i)*this->GetBlockSpacing();
                tmpBlock.SetIndex(startPosition);
                total_num_blocks++;

                if (this->CheckBlockConditions(tmpBlock,blockVar))
                {
                    tmpOutput.push_back(tmpBlock);
                    block_variances.push_back(blockVar);
                }
            }
        }
    }
    else if (NDimensions == 3)
    {
        for (unsigned int k = 0;k < num_blocks[2];++k)
        {
            startPosition[2] = reqStartIndex[2] + (num_start[2] + k)*this->GetBlockSpacing();
            for (unsigned int j = 0;j < num_blocks[1];++j)
            {
                startPosition[1] = reqStartIndex[1] + (num_start[1] + j)*this->GetBlockSpacing();
                for (unsigned int i = 0;i < num_blocks[0];++i)
                {
                    startPosition[0] = reqStartIndex[0] + (num_start[0] + i)*this->GetBlockSpacing();
                    tmpBlock.SetIndex(startPosition);
                    total_num_blocks++;

                    if (this->CheckBlockConditions(tmpBlock,blockVar))
                    {
                        tmpOutput.push_back(tmpBlock);
                        block_variances.push_back(blockVar);
                    }
                }
            }
        }
    }
    else
        throw itk::ExceptionObject(__FILE__, __LINE__,"Images of dimension superior to 3 are not handled...",ITK_LOCATION);
}

template <class PixelType, unsigned int NDimensions>
void BlockMatchingInitializer<PixelType,NDimensions>::ComputeOuterDamFromBlocks()
{
    typedef itk::Image<unsigned char, NDimensions> BlockMaskImageType;
    typedef typename BlockMaskImageType::Pointer BlockMaskImagePointer;

    typedef itk::Image<float, NDimensions> BlockDistanceImageType;

    typename itk::ImageBase<NDimensions>::Pointer refImage = this->GetFirstReferenceImage();

    BlockMaskImagePointer blocksMap = BlockMaskImageType::New();
    blocksMap->Initialize();
    blocksMap->SetRegions (refImage->GetLargestPossibleRegion());
    blocksMap->SetSpacing (refImage->GetSpacing());
    blocksMap->SetOrigin (refImage->GetOrigin());
    blocksMap->SetDirection (refImage->GetDirection());
    blocksMap->Allocate();

    blocksMap->FillBuffer(0);

    IndexType posIndex;
    for (unsigned int i = 0;i < m_Output.size();++i)
    {
        for (unsigned int j = 0;j < NDimensions;++j)
            posIndex[j] = (unsigned int)floor((m_Output[i].GetIndex()[j] + m_Output[i].GetSize()[j] / 2.0) + 0.5);

        blocksMap->SetPixel(posIndex,1);
    }

    typedef itk::DanielssonDistanceMapImageFilter <BlockMaskImageType, BlockDistanceImageType> DistanceMapFilterType;
    typedef typename DistanceMapFilterType::Pointer DistanceMapFilterPointer;

    DistanceMapFilterPointer distMapFilter = DistanceMapFilterType::New();
    distMapFilter->SetInput(blocksMap);
    distMapFilter->SetNumberOfThreads(this->GetNumberOfThreads());
    distMapFilter->InputIsBinaryOn();

    distMapFilter->Update();

    BlockMaskImagePointer damMask = BlockMaskImageType::New();
    damMask->Initialize();
    damMask->SetRegions (refImage->GetLargestPossibleRegion());
    damMask->SetSpacing (refImage->GetSpacing());
    damMask->SetOrigin (refImage->GetOrigin());
    damMask->SetDirection (refImage->GetDirection());
    damMask->Allocate();

    typedef itk::ImageRegionIteratorWithIndex <BlockMaskImageType> MaskIteratorType;
    typedef itk::ImageRegionConstIterator <BlockDistanceImageType> DistMapIteratorType;

    MaskIteratorType damMaskItr(damMask, refImage->GetLargestPossibleRegion());
    DistMapIteratorType distMapItr(distMapFilter->GetOutput(), refImage->GetLargestPossibleRegion());

    IndexType baseIndex, endIndex = refImage->GetLargestPossibleRegion().GetIndex();
    for (unsigned int i = 0;i < NDimensions;++i)
        baseIndex[i] = refImage->GetLargestPossibleRegion().GetSize()[i];

    bool emptyDam = true;
    while (!distMapItr.IsAtEnd())
    {
        double tmpVal = distMapItr.Get();

        // ITK Danielsson distance map is in voxel units
        if ((tmpVal >= m_DamDistance)&&(tmpVal <= m_DamDistance + 1))
        {
            emptyDam = false;
            damMaskItr.Set(1);
            posIndex = damMaskItr.GetIndex();

            for (unsigned int i = 0;i < NDimensions;++i)
            {
                if (baseIndex[i] > posIndex[i])
                    baseIndex[i] = posIndex[i];

                if (endIndex[i] < posIndex[i])
                    endIndex[i] = posIndex[i];
            }
        }
        else
            damMaskItr.Set(0);

        ++distMapItr;
        ++damMaskItr;
    }

    m_DamIndexes.clear();
    if (!emptyDam)
    {
        IndexType tmpIndex;
        // Now select indexes inside bounding box
        if (NDimensions == 2)
        {
            for (unsigned int j = baseIndex[1];j <= endIndex[1];j += m_BlockSpacing)
            {
                tmpIndex[1] = j;
                for (unsigned int i = baseIndex[0];i <= endIndex[0];i += m_BlockSpacing)
                {
                    tmpIndex[0] = i;
                    m_DamIndexes.push_back(tmpIndex);
                }
            }
        }
        else if (NDimensions == 3)
        {
            for (unsigned int k = baseIndex[2];k <= endIndex[2];k += m_BlockSpacing)
            {
                tmpIndex[2] = k;
                for (unsigned int j = baseIndex[1];j <= endIndex[1];j += m_BlockSpacing)
                {
                    tmpIndex[1] = j;
                    for (unsigned int i = baseIndex[0];i <= endIndex[0];i += m_BlockSpacing)
                    {
                        tmpIndex[0] = i;
                        m_DamIndexes.push_back(tmpIndex);
                    }
                }
            }
        }
    }
}

template <class PixelType, unsigned int NDimensions>
bool BlockMatchingInitializer<PixelType,NDimensions>::CheckBlockConditions(ImageRegionType &region, double &blockVariance)
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
        if (!this->CheckTensorVariance(m_ReferenceVectorImages[i],region,tmpVar))
            return false;

        if (tmpVar > blockVariance)
            blockVariance = tmpVar;
    }

    // Reaches here if block respects criterions on all images
    return true;
}

template <class PixelType, unsigned int NDimensions>
bool BlockMatchingInitializer<PixelType,NDimensions>::
CheckTensorVariance(VectorImageType *refImage, ImageRegionType &region, double &blockVariance)
{
    itk::ImageRegionConstIterator <VectorImageType> refItr(refImage,region);
    typedef typename VectorImageType::PixelType VectorType;

    unsigned int nbPts = 0;
    unsigned int vectorSize = refImage->GetNumberOfComponentsPerPixel();

    VectorType meanVal(vectorSize);
    meanVal.Fill(0);
    vnl_matrix <double> blockCovariance(vectorSize,vectorSize);
    blockCovariance.fill(0);

    while (!refItr.IsAtEnd())
    {
        VectorType tmpVal = refItr.Get();
        meanVal += tmpVal;

        for (unsigned int j = 0;j < vectorSize;++j)
            for (unsigned int k = j;k < vectorSize;++k)
                blockCovariance(j,k) += tmpVal[j] * tmpVal[k];

        ++nbPts;

        ++refItr;
    }

    if (nbPts <= 1)
        return false;

    blockVariance = 0;

    for (unsigned int j = 0;j < vectorSize;++j)
        for (unsigned int k = j;k < vectorSize;++k)
        {
            double tmp = (blockCovariance(j,k) - meanVal[j]*meanVal[k]/nbPts)/(nbPts - 1.0);
            if (j == k)
                blockVariance += tmp * tmp;
            else
                blockVariance += 2 * tmp * tmp;
        }

    blockVariance = std::sqrt(blockVariance);

    if (blockVariance > this->GetTensorVarianceThreshold())
        return true;
    else
        return false;
}

template <class PixelType, unsigned int NDimensions>
bool BlockMatchingInitializer<PixelType,NDimensions>::
CheckScalarVariance(ScalarImageType *refImage, ImageRegionType &region, double &blockVariance)
{
    itk::ImageRegionConstIterator <ScalarImageType> refItr(refImage,region);
    blockVariance = 0;
    double meanVal = 0;

    unsigned int nbPts = 0;
    while (!refItr.IsAtEnd())
    {
        double tmpVal = refItr.Get();
        blockVariance += tmpVal*tmpVal;
        meanVal += tmpVal;

        ++nbPts;

        ++refItr;
    }

    if (nbPts <= 1)
        return false;

    blockVariance = (blockVariance - meanVal * meanVal / nbPts) / (nbPts - 1.0);

    if (blockVariance > this->GetScalarVarianceThreshold())
        return true;
    else
        return false;
}

}// end of namespace anima
