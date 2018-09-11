#pragma once
#include "animaBaseBlockMatcher.h"

#include <animaBobyqaOptimizer.h>
#include <animaVoxelExhaustiveOptimizer.h>
#include <animaBlockMatchInitializer.h>

namespace anima
{

template <typename TInputImageType>
BaseBlockMatcher <TInputImageType>
::BaseBlockMatcher()
{
    m_ForceComputeBlocks = false;
    m_NumberOfThreads = itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads();

    m_BlockPercentageKept = 0.8;
    m_BlockSize = 5;
    m_BlockSpacing = 2;

    m_BlockVarianceThreshold = 25;

    m_SearchRadius = 2;
    m_FinalRadius = 0.001;
    m_OptimizerMaximumIterations = 100;
    m_StepSize = 1.0;

    m_OptimizerType = Bobyqa;
    m_Verbose = true;

    m_HighestProcessedBlock = 0;
}

template <typename TInputImageType>
void
BaseBlockMatcher <TInputImageType>
::InitializeBlocks()
{
    // Init blocks on reference image
    typedef typename TInputImageType::IOPixelType InputPixelType;
    typedef typename anima::BlockMatchingInitializer<InputPixelType,TInputImageType::ImageDimension> InitializerType;
    typedef typename InitializerType::Pointer InitializerPointer;

    InitializerPointer initPtr = InitializerType::New();
    initPtr->AddReferenceImage(m_ReferenceImage);

    if (m_NumberOfThreads != 0)
        initPtr->SetNumberOfThreads(m_NumberOfThreads);

    initPtr->SetPercentageKept(m_BlockPercentageKept);
    initPtr->SetBlockSize(m_BlockSize);
    initPtr->SetBlockSpacing(m_BlockSpacing);
    initPtr->SetScalarVarianceThreshold(m_BlockVarianceThreshold);
    initPtr->SetOrientedModelVarianceThreshold(m_BlockVarianceThreshold);
    initPtr->AddGenerationMask(m_BlockGenerationMask);

    initPtr->SetRequestedRegion(m_ReferenceImage->GetLargestPossibleRegion());

    m_BlockRegions = initPtr->GetOutput();
    m_BlockPositions = initPtr->GetOutputPositions();

    if (m_Verbose)
        std::cout << "Generated " << m_BlockRegions.size() << " blocks..." << std::endl;

    m_BlockTransformPointers.resize(m_BlockRegions.size());
    m_BlockWeights.resize(m_BlockRegions.size());
    for (unsigned int i = 0;i < m_BlockRegions.size();++i)
        m_BlockTransformPointers[i] = this->GetNewBlockTransform(m_BlockPositions[i]);
}

template <typename TInputImageType>
typename BaseBlockMatcher <TInputImageType>::OptimizerPointer
BaseBlockMatcher <TInputImageType>
::SetupOptimizer()
{
    OptimizerPointer optimizer;
    switch (m_OptimizerType)
    {
        case Bobyqa:
        {
            typedef anima::BobyqaOptimizer LocalOptimizerType;
            optimizer = LocalOptimizerType::New();
            LocalOptimizerType *tmpOpt = (LocalOptimizerType *)optimizer.GetPointer();
            tmpOpt->SetRhoBegin(m_SearchRadius);
            tmpOpt->SetRhoEnd(m_FinalRadius);

            tmpOpt->SetNumberSamplingPoints(m_BlockTransformPointers[0]->GetNumberOfParameters() + 2);
            tmpOpt->SetMaximumIteration(m_OptimizerMaximumIterations);
            tmpOpt->SetMaximize(this->GetMaximizedMetric());

            break;
        }

        case Exhaustive:
        default:
        {
            typedef anima::VoxelExhaustiveOptimizer LocalOptimizerType;
            optimizer = LocalOptimizerType::New();

            LocalOptimizerType *tmpOpt = (LocalOptimizerType *)optimizer.GetPointer();
            LocalOptimizerType::StepsType steps(m_BlockTransformPointers[0]->GetNumberOfParameters());

            typename InputImageType::SpacingType fixedSpacing = m_ReferenceImage->GetSpacing();
            typename InputImageType::DirectionType fixedDirection = m_ReferenceImage->GetDirection();

            vnl_matrix <double> geometry(InputImageType::ImageDimension,InputImageType::ImageDimension);
            for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
            {
                double tmpVoxSize = fixedSpacing[i];

                for (unsigned int j = 0;j < InputImageType::ImageDimension;++j)
                    geometry(j,i) = tmpVoxSize*fixedDirection(j,i);
            }

            tmpOpt->SetGeometry(geometry);

            LocalOptimizerType::ScalesType tmpScales(InputImageType::ImageDimension);

            for (unsigned i = 0; i < steps.Size(); ++i)
            {
                steps[i] = m_SearchRadius;
                tmpScales[i] = m_StepSize;
            }

            tmpOpt->SetNumberOfSteps(steps);
            tmpOpt->SetScales(tmpScales);
            tmpOpt->SetMaximize(this->GetMaximizedMetric());
            break;
        }
    }

    this->TransformDependantOptimizerSetup(optimizer);

    return optimizer;
}

template <typename TInputImageType>
void
BaseBlockMatcher <TInputImageType>
::Update()
{
    // Generate blocks if needed on reference image
    if ((m_ForceComputeBlocks) || (m_BlockTransformPointers.size() == 0))
        this->InitializeBlocks();

    m_HighestProcessedBlock = 0;
    itk::MultiThreaderBase::Pointer threadWorker = itk::MultiThreaderBase::New();
    ThreadedMatchData *tmpStr = new ThreadedMatchData;
    tmpStr->BlockMatch = this;

    threadWorker->SetNumberOfWorkUnits(m_NumberOfThreads);
    threadWorker->SetSingleMethod(this->ThreadedMatching,tmpStr);
    threadWorker->SingleMethodExecute();

    delete tmpStr;
}

template <typename TInputImageType>
ITK_THREAD_RETURN_TYPE
BaseBlockMatcher <TInputImageType>
::ThreadedMatching(void *arg)
{
    itk::MultiThreaderBase::WorkUnitInfo *threadArgs = (itk::MultiThreaderBase::WorkUnitInfo *)arg;
    ThreadedMatchData* data = (ThreadedMatchData *)threadArgs->UserData;

    data->BlockMatch->ProcessBlockMatch();
    return NULL;
}

template <typename TInputImageType>
void
BaseBlockMatcher <TInputImageType>
::ProcessBlockMatch()
{
    bool continueLoop = true;
    unsigned int highestToleratedBlockIndex = m_BlockRegions.size();

    unsigned int stepData = std::min((int)m_BlockRegions.size(),100);
    if (stepData == 0)
        stepData = 1;

    while (continueLoop)
    {
        m_LockHighestProcessedBlock.Lock();

        if (m_HighestProcessedBlock >= highestToleratedBlockIndex)
        {
            m_LockHighestProcessedBlock.Unlock();
            continueLoop = false;
            continue;
        }

        unsigned int startPoint = m_HighestProcessedBlock;
        unsigned int endPoint = m_HighestProcessedBlock + stepData;
        if (endPoint > highestToleratedBlockIndex)
            endPoint = highestToleratedBlockIndex;

        m_HighestProcessedBlock = endPoint;

        m_LockHighestProcessedBlock.Unlock();

        this->BlockMatch(startPoint,endPoint);
    }
}

template <typename TInputImageType>
void
BaseBlockMatcher <TInputImageType>
::BlockMatch(unsigned int startIndex, unsigned int endIndex)
{
    MetricPointer metric = this->SetupMetric();
    OptimizerPointer optimizer = this->SetupOptimizer();

    // Loop over the desired blocks
    for (unsigned int block = startIndex;block < endIndex;++block)
    {
        this->BlockMatchingSetup(metric, block);
        optimizer->SetCostFunction(metric);
        optimizer->SetInitialPosition(m_BlockTransformPointers[block]->GetParameters());

        try
        {
            optimizer->StartOptimization();
        }
        catch (itk::ExceptionObject & err)
        {
            m_BlockWeights[block] = 0;
            continue;
        }

        m_BlockTransformPointers[block]->SetParameters(optimizer->GetCurrentPosition());

        double val = optimizer->GetValue(optimizer->GetCurrentPosition());
        m_BlockWeights[block] = this->ComputeBlockWeight(val,block);
    }
}

} // end namespace anima
