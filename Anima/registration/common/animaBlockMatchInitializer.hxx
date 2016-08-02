#pragma once
#include "animaBlockMatchInitializer.h"

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageFileWriter.h>

#include <itkExpNegativeImageFilter.h>
#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>

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
typename BlockMatchingInitializer<PixelType,NDimensions>::WeightImagePointer &
BlockMatchingInitializer<PixelType,NDimensions>
::GetBlockDamWeights()
{
    this->Update();

    return m_BlockDamWeights;
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
::SetTensorVarianceThreshold(double val)
{
    if (val != m_TensorVarianceThreshold)
    {
        m_TensorVarianceThreshold = val;
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
::SetComputeOuterDam(bool val)
{
    if (val != m_ComputeOuterDam)
    {
        m_ComputeOuterDam = val;
        m_UpToDate = false;
    }
}

template <class PixelType, unsigned int NDimensions>
void BlockMatchingInitializer<PixelType,NDimensions>
::SetDamDistance(double val)
{
    if (val != m_DamDistance)
    {
        m_DamDistance = val;
        m_UpToDate = false;
    }
}

template <class PixelType, unsigned int NDimensions>
void BlockMatchingInitializer<PixelType,NDimensions>
::Update()
{
    if (m_UpToDate)
        return;

    if (m_ComputeOuterDam)
    {
        m_BlockDamWeights = WeightImageType::New();

        m_BlockDamWeights->Initialize();
        m_BlockDamWeights->SetRegions (this->GetFirstReferenceImage()->GetLargestPossibleRegion());
        m_BlockDamWeights->SetSpacing (this->GetFirstReferenceImage()->GetSpacing());
        m_BlockDamWeights->SetOrigin (this->GetFirstReferenceImage()->GetOrigin());
        m_BlockDamWeights->SetDirection (this->GetFirstReferenceImage()->GetDirection());
        m_BlockDamWeights->Allocate();

        m_BlockDamWeights->FillBuffer(0);
    }

    itk::MultiThreader::Pointer threaderBlockGenerator = itk::MultiThreader::New();

    BlockGeneratorThreadStruct *tmpStr = new BlockGeneratorThreadStruct;

    std::vector <unsigned int> totalNbBlocks(NDimensions);
    tmpStr->blockStartOffsets.resize(NDimensions);

    for (unsigned int i = 0;i < NDimensions;++i)
    {
        totalNbBlocks[i] = floor((float)(this->GetRequestedRegion().GetSize()[i] / this->GetBlockSpacing()) + 0.5);
        unsigned int spaceRequired = this->GetRequestedRegion().GetSize()[i] - (totalNbBlocks[i] - 1) * this->GetBlockSpacing() - 1;

        tmpStr->blockStartOffsets[i] = floor(spaceRequired / 2.0);
    }

    unsigned int nb_blocks_per_thread = (unsigned int) floor((float)(totalNbBlocks[NDimensions-1] / this->GetNumberOfThreads()));
    if (nb_blocks_per_thread < 1)
    {
        nb_blocks_per_thread = 1;
        this->SetNumberOfThreads(totalNbBlocks[NDimensions-1]);
    }

    tmpStr->startBlocks.resize(this->GetNumberOfThreads());
    tmpStr->nb_blocks.resize(this->GetNumberOfThreads());

    unsigned int currentCount = 0;
    for (unsigned int i = 0;i < this->GetNumberOfThreads();++i)
    {
        std::vector <unsigned int> startBlock(NDimensions,0);
        startBlock[NDimensions-1] = currentCount;
        tmpStr->startBlocks[i] = startBlock;

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
        tmpStr->nb_blocks[i] = numBlockVec;

        currentCount += numBlocksForThread;
    }

    tmpStr->Filter = this;
    tmpStr->tmpOutput.resize(this->GetNumberOfThreads());
    tmpStr->totalNumberOfBlocks.resize(this->GetNumberOfThreads());
    tmpStr->blocks_positions.resize(this->GetNumberOfThreads());
    tmpStr->blocks_variances.resize(this->GetNumberOfThreads());

    for (unsigned int i = 0;i < this->GetNumberOfThreads();++i)
    {
        tmpStr->tmpOutput[i].clear();
        tmpStr->totalNumberOfBlocks[i] = 0;
        tmpStr->blocks_positions[i].clear();
        tmpStr->blocks_variances[i].clear();
    }

    threaderBlockGenerator->SetNumberOfThreads(this->GetNumberOfThreads());
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

    if (m_ComputeOuterDam)
        this->ComputeOuterDamFromBlocks();

    m_UpToDate = true;
}

template <class PixelType, unsigned int NDimensions>
ITK_THREAD_RETURN_TYPE
BlockMatchingInitializer<PixelType,NDimensions>
::ThreadBlockGenerator(void *arg)
{
    itk::MultiThreader::ThreadInfoStruct *threadArgs = (itk::MultiThreader::ThreadInfoStruct *)arg;

    unsigned int nbThread = threadArgs->ThreadID;

    BlockGeneratorThreadStruct *tmpStr = (BlockGeneratorThreadStruct *)threadArgs->UserData;

    tmpStr->Filter->RegionBlockGenerator(tmpStr->startBlocks[nbThread], tmpStr->blockStartOffsets, tmpStr->nb_blocks[nbThread],
                                         tmpStr->tmpOutput[nbThread], tmpStr->blocks_positions[nbThread],
                                         tmpStr->blocks_variances[nbThread], tmpStr->totalNumberOfBlocks[nbThread]);

    return NULL;
}

template <class PixelType, unsigned int NDimensions>
void BlockMatchingInitializer<PixelType,NDimensions>
::RegionBlockGenerator(std::vector <unsigned int> &num_start, std::vector <unsigned int> &block_start_offsets,
                       std::vector <unsigned int> &num_blocks, std::vector <ImageRegionType> &tmpOutput,
                       std::vector<PointType> &block_origins, std::vector <double> &block_variances,
                       unsigned int &total_num_blocks)
{
    typename ImageRegionType::IndexType startPosition, blockPosition;
    typename ImageRegionType::SizeType blockSize;
    PointType blockOrigin;
    ImageRegionType largestRegion = this->GetFirstReferenceImage()->GetLargestPossibleRegion();

    ImageRegionType tmpBlock;
    int indexPos;
    double blockVar = 0;
    total_num_blocks = 0;
    block_variances.clear();
    block_origins.clear();

    unsigned int block_half_size = floor ((this->GetBlockSize() - 1) / 2.0);

    bool continueLoop = true;
    std::vector <unsigned int> positionCounter(NDimensions,0);

    while (continueLoop)
    {
        for (unsigned int i = 0;i < NDimensions;++i)
        {
            blockSize[i] = this->GetBlockSize();
            indexPos = block_start_offsets[i] + (num_start[i] + positionCounter[i])*this->GetBlockSpacing() - block_half_size;
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
        total_num_blocks++;

        if (this->CheckBlockConditions(tmpBlock,blockVar))
        {
            tmpOutput.push_back(tmpBlock);
            block_variances.push_back(blockVar);
            for (unsigned int i = 0;i < NDimensions;++i)
                blockPosition[i] = block_start_offsets[i] + (num_start[i] + positionCounter[i])*this->GetBlockSpacing();
            this->GetFirstReferenceImage()->TransformIndexToPhysicalPoint(blockPosition,blockOrigin);
            block_origins.push_back(blockOrigin);
        }

        continueLoop = this->ProgressCounter(positionCounter,num_blocks);
    }
}

template <class PixelType, unsigned int NDimensions>
void BlockMatchingInitializer<PixelType,NDimensions>
::ComputeOuterDamFromBlocks()
{
    IndexType posIndex;
    for (unsigned int i = 0;i < m_Output.size();++i)
    {
        for (unsigned int j = 0;j < NDimensions;++j)
            posIndex[j] = (unsigned int)std::round(m_Output[i].GetIndex()[j] + (m_Output[i].GetSize()[j] - 1) / 2.0);

        m_BlockDamWeights->SetPixel(posIndex,1);
    }

    typedef itk::BinaryBallStructuringElement <unsigned short, NDimensions> BallElementType;
    typedef itk::GrayscaleDilateImageFilter <WeightImageType,WeightImageType,BallElementType> DilateFilterType;
    typename DilateFilterType::Pointer dilateFilter = DilateFilterType::New();

    dilateFilter->SetInput(m_BlockDamWeights);

    BallElementType dilateBall;
    // Rule of thumb for dilation of block mask
    unsigned int radius = std::floor(m_BlockSpacing + (m_BlockSize - 1.0) / 2.0);
    dilateBall.SetRadius(radius);
    dilateBall.CreateStructuringElement();

    dilateFilter->SetKernel(dilateBall);
    dilateFilter->SetNumberOfThreads(this->GetNumberOfThreads());

    dilateFilter->Update();

    typedef itk::DanielssonDistanceMapImageFilter <WeightImageType, WeightImageType> DistanceMapFilterType;
    typedef typename DistanceMapFilterType::Pointer DistanceMapFilterPointer;

    DistanceMapFilterPointer distMapFilter = DistanceMapFilterType::New();
    distMapFilter->SetInput(dilateFilter->GetOutput());
    distMapFilter->SetNumberOfThreads(this->GetNumberOfThreads());
    distMapFilter->InputIsBinaryOn();
    distMapFilter->SetSquaredDistance(true);
    distMapFilter->SetUseImageSpacing(true);

    distMapFilter->Update();

    typedef itk::ExpNegativeImageFilter <WeightImageType,WeightImageType> WeightFilterType;
    typename WeightFilterType::Pointer weightFilter = WeightFilterType::New();

    weightFilter->SetInput(distMapFilter->GetOutput());
    weightFilter->SetFactor(1.0 / (m_DamDistance * m_DamDistance));
    weightFilter->SetNumberOfThreads(this->GetNumberOfThreads());

    weightFilter->Update();

    m_BlockDamWeights = weightFilter->GetOutput();
    m_BlockDamWeights->DisconnectPipeline();
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
::CheckBlockConditions(ImageRegionType &region, double &blockVariance)
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
bool BlockMatchingInitializer<PixelType,NDimensions>
::CheckTensorVariance(VectorImageType *refImage, ImageRegionType &region, double &blockVariance)
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
