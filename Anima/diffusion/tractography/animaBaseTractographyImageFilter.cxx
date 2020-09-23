#include "animaBaseTractographyImageFilter.h"

#include <itkImageRegionIteratorWithIndex.h>
#include <itkMultiThreaderBase.h>
#include <itkProgressReporter.h>

#include <animaVectorOperations.h>

#include <vtkPointData.h>
#include <vtkCellData.h>

namespace anima
{

BaseTractographyImageFilter::BaseTractographyImageFilter()
{
    m_PointsToProcess.clear();

    m_NumberOfFibersPerPixel = 1;
    
    m_StepProgression = 1;
    m_MaxFiberAngle = M_PI / 2.0;

    m_ComputeLocalColors = true;
    m_HighestProcessedSeed = 0;
    m_ProgressReport = ITK_NULLPTR;
}

BaseTractographyImageFilter::~BaseTractographyImageFilter()
{
    if (m_ProgressReport)
        delete m_ProgressReport;
}

void BaseTractographyImageFilter::Update()
{
    this->PrepareTractography();
    m_Output = vtkPolyData::New();
    
    if (m_ProgressReport)
        delete m_ProgressReport;

    unsigned int stepData = std::min((int)m_PointsToProcess.size(),100);
    if (stepData == 0)
        stepData = 1;

    unsigned int numSteps = std::floor(m_PointsToProcess.size() / (double)stepData);
    if (m_PointsToProcess.size() % stepData != 0)
        numSteps++;

    m_ProgressReport = new itk::ProgressReporter(this,0,numSteps);

    std::vector < FiberType > resultFibers;
    
    trackerArguments tmpStr;
    tmpStr.trackerPtr = this;
    for (unsigned int i = 0;i < this->GetNumberOfWorkUnits();++i)
        tmpStr.resultFibersFromThreads.push_back(resultFibers);

    this->GetMultiThreader()->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    this->GetMultiThreader()->SetSingleMethod(this->ThreadTracker,&tmpStr);
    this->GetMultiThreader()->SingleMethodExecute();
    
    for (unsigned int j = 0;j < this->GetNumberOfWorkUnits();++j)
    {
        resultFibers.insert(resultFibers.end(),tmpStr.resultFibersFromThreads[j].begin(),
                            tmpStr.resultFibersFromThreads[j].end());
    }
    
    std::cout << "\nTracked a total of " << resultFibers.size() << " fibers" << std::endl;
    resultFibers = this->FilterOutputFibers(resultFibers);
    std::cout << "Kept " << resultFibers.size() << " fibers after filtering" << std::endl;
    this->createVTKOutput(resultFibers);
}

void BaseTractographyImageFilter::PrepareTractography()
{
    typedef itk::ImageRegionIteratorWithIndex <MaskImageType> MaskImageIteratorType;
    
    MaskImageIteratorType seedItr(m_SeedingImage, m_InputImage->GetLargestPossibleRegion());
    m_PointsToProcess.clear();
    
    IndexType tmpIndex;
    PointType tmpPoint;
    ContinuousIndexType realIndex;
    
    m_FilteringValues.clear();
    double startN = -0.5 + 1.0 / (2.0 * m_NumberOfFibersPerPixel);
    double stepN = 1.0 / m_NumberOfFibersPerPixel;
    
    bool is2d = (m_InputImage->GetLargestPossibleRegion().GetSize()[2] <= 1);
    FiberType tmpFiber(1);
    
    while (!seedItr.IsAtEnd())
    {
        if (seedItr.Get() == 0)
        {
            ++seedItr;
            continue;
        }
        
        tmpIndex = seedItr.GetIndex();
        
        if (is2d)
        {
            realIndex[2] = tmpIndex[2];
            for (unsigned int j = 0;j < m_NumberOfFibersPerPixel;++j)
            {
                realIndex[1] = tmpIndex[1] + startN + j * stepN;
                for (unsigned int i = 0;i < m_NumberOfFibersPerPixel;++i)
                {
                    realIndex[0] = tmpIndex[0] + startN + i * stepN;
                    m_SeedingImage->TransformContinuousIndexToPhysicalPoint(realIndex,tmpPoint);
                    tmpFiber[0] = tmpPoint;
                    m_PointsToProcess.push_back(tmpFiber);
                }
            }
        }
        else
        {
            for (unsigned int k = 0;k < m_NumberOfFibersPerPixel;++k)
            {
                realIndex[2] = tmpIndex[2] + startN + k * stepN;
                for (unsigned int j = 0;j < m_NumberOfFibersPerPixel;++j)
                {
                    realIndex[1] = tmpIndex[1] + startN + j * stepN;
                    for (unsigned int i = 0;i < m_NumberOfFibersPerPixel;++i)
                    {
                        realIndex[0] = tmpIndex[0] + startN + i * stepN;
                        m_SeedingImage->TransformContinuousIndexToPhysicalPoint(realIndex,tmpPoint);
                        tmpFiber[0] = tmpPoint;
                        m_PointsToProcess.push_back(tmpFiber);
                    }
                }
            }
        }
        
        ++seedItr;
    }

    if (m_FilteringImage.IsNotNull())
    {
        MaskImageIteratorType filterItr(m_FilteringImage, m_InputImage->GetLargestPossibleRegion());

        while (!filterItr.IsAtEnd())
        {
            if (filterItr.Get() == 0)
            {
                ++filterItr;
                continue;
            }

            bool isAlreadyIn = false;
            for (unsigned int i = 0;i < m_FilteringValues.size();++i)
            {
                if (m_FilteringValues[i] == filterItr.Get())
                {
                    isAlreadyIn = true;
                    break;
                }
            }

            if (!isAlreadyIn)
                m_FilteringValues.push_back(filterItr.Get());

            ++filterItr;
        }
    }
    
    std::cout << "Generated " << m_PointsToProcess.size() << " seed points from ROI mask" << std::endl;
}

ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION BaseTractographyImageFilter::ThreadTracker(void *arg)
{
    itk::MultiThreaderBase::WorkUnitInfo *threadArgs = (itk::MultiThreaderBase::WorkUnitInfo *)arg;
    unsigned int nbThread = threadArgs->WorkUnitID;
    
    trackerArguments *tmpArg = (trackerArguments *)threadArgs->UserData;
    tmpArg->trackerPtr->ThreadTrack(nbThread,tmpArg->resultFibersFromThreads[nbThread]);
    
    return ITK_THREAD_RETURN_DEFAULT_VALUE;
}

void BaseTractographyImageFilter::ThreadTrack(unsigned int numThread, std::vector <FiberType> &resultFibers)
{
    bool continueLoop = true;
    unsigned int highestToleratedSeedIndex = m_PointsToProcess.size();

    unsigned int stepData = std::min((int)m_PointsToProcess.size(),100);
    if (stepData == 0)
        stepData = 1;

    while (continueLoop)
    {
        m_LockHighestProcessedSeed.lock();

        if (m_HighestProcessedSeed >= (int)highestToleratedSeedIndex)
        {
            m_LockHighestProcessedSeed.unlock();
            continueLoop = false;
            continue;
        }

        unsigned int startPoint = m_HighestProcessedSeed;
        unsigned int endPoint = m_HighestProcessedSeed + stepData;
        if (endPoint > highestToleratedSeedIndex)
            endPoint = highestToleratedSeedIndex;

        m_HighestProcessedSeed = endPoint;

        m_LockHighestProcessedSeed.unlock();

        this->ThreadedTrackComputer(numThread,resultFibers,startPoint,endPoint);

        m_LockHighestProcessedSeed.lock();
        m_ProgressReport->CompletedPixel();
        m_LockHighestProcessedSeed.unlock();
    }
}

void BaseTractographyImageFilter::ThreadedTrackComputer(unsigned int numThread, std::vector <FiberType> &resultFibers,
                                                        unsigned int startSeedIndex, unsigned int endSeedIndex)
{    
    std::vector <PointType> initialDirections;
    VectorType modelValue;
    bool is2d = (m_SeedingImage->GetLargestPossibleRegion().GetSize()[2] <= 1);
    FiberType tmpFiber;
    ContinuousIndexType curIndex;
    IndexType curNearestIndex;

    for (unsigned int i = startSeedIndex;i < endSeedIndex;++i)
    {
        bool treatPoint = true;

        m_SeedingImage->TransformPhysicalPointToContinuousIndex(m_PointsToProcess[i][0],curIndex);
        m_SeedingImage->TransformPhysicalPointToIndex(m_PointsToProcess[i][0],curNearestIndex);
        modelValue.SetSize(1);

        if (!m_CutMaskImage.IsNull())
        {
            if (!this->CheckIndexInImageBounds(curNearestIndex,m_CutMaskImage))
                treatPoint = false;
            else if (m_CutMaskImage->GetPixel(curNearestIndex) != 0)
                treatPoint = false;
        }

        if (!m_ForbiddenMaskImage.IsNull())
        {
            if (!this->CheckIndexInImageBounds(curNearestIndex,m_ForbiddenMaskImage))
                treatPoint = false;
            else if (m_ForbiddenMaskImage->GetPixel(curNearestIndex) != 0)
                treatPoint = false;
        }

        if (!treatPoint)
            continue;

        this->GetModelValue(curIndex,modelValue);
        if (isZero(modelValue))
            treatPoint = false;

        if (!this->CheckModelCompatibility(modelValue,numThread))
            treatPoint = false;

        initialDirections = this->GetModelPrincipalDirections(modelValue, is2d, numThread);

        unsigned int numDirections = initialDirections.size();
        for (unsigned int j = 0;j < numDirections;++j)
        {
            tmpFiber = m_PointsToProcess[i];
            this->ComputeFiber(tmpFiber,initialDirections[j],numThread);
        
            if (tmpFiber.size() > m_MinLengthFiber / m_StepProgression)
                resultFibers.push_back(tmpFiber);
        }
    }
}

std::vector < std::vector <BaseTractographyImageFilter::PointType> >
BaseTractographyImageFilter::FilterOutputFibers(std::vector < FiberType > &fibers)
{
    std::vector < FiberType > resVal;
    
    if (m_FilteringValues.size() > 1)
    {
        std::vector <bool> touchingLabels(m_FilteringValues.size());
        IndexType tmpIndex;
        PointType tmpPoint;
        for (unsigned int i = 0;i < fibers.size();++i)
        {
            std::fill(touchingLabels.begin(),touchingLabels.end(),false);
            for (unsigned int j = 0;j < fibers[i].size();++j)
            {
                tmpPoint = fibers[i][j];
                m_FilteringImage->TransformPhysicalPointToIndex(tmpPoint,tmpIndex);
                
                unsigned int maskValue = m_FilteringImage->GetPixel(tmpIndex);
                if (maskValue != 0)
                {
                    for (unsigned int k = 0;k < m_FilteringValues.size();++k)
                    {
                        if (maskValue == m_FilteringValues[k])
                        {
                            touchingLabels[k] = true;
                            break;
                        }
                    }
                }
            }
            
            bool keepFiber = true;
            for (unsigned int k = 0;k < touchingLabels.size();++k)
            {
                if (!touchingLabels[k])
                {
                    keepFiber = false;
                    break;
                }
            }
            
            if (keepFiber)
                resVal.push_back(fibers[i]);
        }
    }
    else
        resVal = fibers;
    
    return resVal;
}

void BaseTractographyImageFilter::createVTKOutput(std::vector < FiberType > &filteredFibers)
{
    m_Output = vtkPolyData::New();
    m_Output->Initialize();
    m_Output->Allocate();
    
    vtkSmartPointer <vtkPoints> myPoints = vtkPoints::New();    
    for (unsigned int i = 0;i < filteredFibers.size();++i)
    {
        unsigned int npts = filteredFibers[i].size();
        vtkIdType* ids = new vtkIdType[npts];
        
        for (unsigned int j = 0;j < npts;++j)
            ids[j] = myPoints->InsertNextPoint(filteredFibers[i][j][0],filteredFibers[i][j][1],filteredFibers[i][j][2]);
        
        m_Output->InsertNextCell (VTK_POLY_LINE, npts, ids);
        delete[] ids;
    }
    
    m_Output->SetPoints (myPoints);
    if (m_ComputeLocalColors)
    {
        std::cout << "Computing local colors and microstructure maps" << std::endl;
        this->ComputeAdditionalScalarMaps();
    }
}

BaseTractographyImageFilter::FiberVectorType
BaseTractographyImageFilter::ComputeFiber(BaseTractographyImageFilter::FiberType &fiber,
                                          PointType &initialDirection,
                                          itk::ThreadIdType threadId)
{
    FiberVectorType resVal;
    bool is2d = m_InputImage->GetLargestPossibleRegion().GetSize()[2] == 1;
    
    PointType curPoint = fiber[0];
    PointType oldDir, newDir;
    oldDir = initialDirection;
    newDir = initialDirection;

    ContinuousIndexType curIndex;
    IndexType curNearestIndex;
    bool continueLoop = true;
    VectorType modelValue;
    modelValue.SetSize(1);

    this->GetModelValue(curIndex,modelValue);
    // Go the forward way starting along initial direction provided
    while (continueLoop)
    {
        if (fiber.size() > 1)
        {
            for (unsigned int i = 0;i < 3;++i)
                oldDir[i] = fiber[fiber.size() - 1][i] - fiber[fiber.size() - 2][i];

            newDir = this->GetNextDirection(oldDir,modelValue,is2d,threadId);
        }

        if (anima::ComputeOrientationAngle(oldDir, newDir) > m_MaxFiberAngle)
        {
            continueLoop = false;
            break;
        }

        this->ComputeNewFiberPoint(fiber[fiber.size() - 1],newDir,curPoint,threadId);

        m_SeedingImage->TransformPhysicalPointToContinuousIndex(curPoint,curIndex);
        m_SeedingImage->TransformPhysicalPointToIndex(curPoint,curNearestIndex);

        if (!m_CutMaskImage.IsNull())
        {
            if (!this->CheckIndexInImageBounds(curNearestIndex,m_CutMaskImage))
                continueLoop = false;
            else if (m_CutMaskImage->GetPixel(curNearestIndex) != 0)
                continueLoop = false;
        }

        if (!m_ForbiddenMaskImage.IsNull())
        {
            if (!this->CheckIndexInImageBounds(curNearestIndex,m_ForbiddenMaskImage))
            {
                fiber.clear();
                return resVal;
            }
            if (m_ForbiddenMaskImage->GetPixel(curNearestIndex) != 0)
            {
                fiber.clear();
                return resVal;
            }
        }

        if (!this->CheckIndexInImageBounds(curIndex))
            continueLoop = false;
        else
        {
            this->GetModelValue(curIndex,modelValue);

            if (isZero(modelValue))
                continueLoop = false;
        }

        if (!this->CheckModelCompatibility(modelValue,threadId))
            continueLoop = false;

        // Add new point to fiber
        if (continueLoop)
            fiber.push_back(curPoint);

        if (fiber.size() > m_MaxLengthFiber / m_StepProgression)
        {
            fiber.clear();
            return resVal;
        }
    }

    curPoint = fiber[0];
    oldDir = initialDirection;
    for (unsigned int i = 0;i < 3;++i)
        oldDir[i] *= -1.0;

    newDir = oldDir;
    continueLoop = true;
    m_SeedingImage->TransformPhysicalPointToContinuousIndex(curPoint,curIndex);
    m_SeedingImage->TransformPhysicalPointToIndex(curPoint,curNearestIndex);

    modelValue.SetSize(1);
    this->GetModelValue(curIndex,modelValue);

    // Now go the backwards way to complete the fiber
    while (continueLoop)
    {
        if (fiber.size() > 1)
        {
            for (unsigned int i = 0;i < 3;++i)
                oldDir[i] = fiber[0][i] - fiber[1][i];

            newDir = this->GetNextDirection(oldDir,modelValue,is2d,threadId);
        }

        if (anima::ComputeOrientationAngle(oldDir, newDir) > m_MaxFiberAngle)
        {
            continueLoop = false;
            break;
        }

        this->ComputeNewFiberPoint(fiber[0],newDir,curPoint,threadId);

        m_SeedingImage->TransformPhysicalPointToContinuousIndex(curPoint,curIndex);
        m_SeedingImage->TransformPhysicalPointToIndex(curPoint,curNearestIndex);

        if (!m_CutMaskImage.IsNull())
        {
            if (!this->CheckIndexInImageBounds(curNearestIndex,m_CutMaskImage))
                continueLoop = false;
            else if (m_CutMaskImage->GetPixel(curNearestIndex) != 0)
                continueLoop = false;
        }

        if (!m_ForbiddenMaskImage.IsNull())
        {
            if (!this->CheckIndexInImageBounds(curNearestIndex,m_ForbiddenMaskImage))
            {
                fiber.clear();
                return resVal;
            }
            if (m_ForbiddenMaskImage->GetPixel(curNearestIndex) != 0)
            {
                fiber.clear();
                return resVal;
            }
        }

        if (!this->CheckIndexInImageBounds(curIndex))
            continueLoop = false;
        else
        {
            this->GetModelValue(curIndex,modelValue);

            if (isZero(modelValue))
                continueLoop = false;
        }

        if (!this->CheckModelCompatibility(modelValue,threadId))
            continueLoop = false;

        // Now add point to fiber
        if (continueLoop)
            fiber.insert(fiber.begin(),1,curPoint);

        if (fiber.size() > m_MaxLengthFiber / m_StepProgression)
        {
            fiber.clear();
            return resVal;
        }
    }
    
    return resVal;
}

void BaseTractographyImageFilter::ComputeNewFiberPoint(PointType &oldPoint, PointType &newDirection, PointType &newPoint,
                                                       itk::ThreadIdType threadId)
{
    PointType k1, k2, k3, k4;
    PointType kPosition;
    ContinuousIndexType curIndex;
    VectorType modelValue;

    bool is2d = (m_InputImage->GetLargestPossibleRegion().GetSize()[2] <= 1);

    for (unsigned int i = 0;i < 3;++i)
    {
        k1[i] = newDirection[i];
        kPosition[i] = oldPoint[i] + m_StepProgression * k1[i] / 2.0;
    }

    m_SeedingImage->TransformPhysicalPointToContinuousIndex(kPosition,curIndex);
    if (!this->CheckIndexInImageBounds(curIndex))
    {
        for (unsigned int i = 0;i < 3;++i)
            newPoint[i] = kPosition[i];

        return;
    }

    this->GetModelValue(curIndex,modelValue);
    if (isZero(modelValue))
    {
        for (unsigned int i = 0;i < 3;++i)
            newPoint[i] = kPosition[i];

        return;
    }

    k2 = this->GetNextDirection(k1,modelValue,is2d,threadId);
    for (unsigned int i = 0;i < 3;++i)
        kPosition[i] = oldPoint[i] + m_StepProgression * k2[i] / 2.0;

    m_SeedingImage->TransformPhysicalPointToContinuousIndex(kPosition,curIndex);
    if (!this->CheckIndexInImageBounds(curIndex))
    {
        for (unsigned int i = 0;i < 3;++i)
            newPoint[i] = oldPoint[i] + (k1[i] + 2.0 * k2[i]) / 3.0;

        return;
    }

    this->GetModelValue(curIndex,modelValue);
    if (isZero(modelValue))
    {
        for (unsigned int i = 0;i < 3;++i)
            newPoint[i] = oldPoint[i] + (k1[i] + 2.0 * k2[i]) / 3.0;

        return;
    }

    k3 = this->GetNextDirection(k2,modelValue,is2d,threadId);
    for (unsigned int i = 0;i < 3;++i)
        kPosition[i] = oldPoint[i] + m_StepProgression * k3[i];

    m_SeedingImage->TransformPhysicalPointToContinuousIndex(kPosition,curIndex);
    if (!this->CheckIndexInImageBounds(curIndex))
    {
        for (unsigned int i = 0;i < 3;++i)
            newPoint[i] = oldPoint[i] + (k1[i] + 2.0 * k2[i] + 2.0 * k3[i]) / 5.0;

        return;
    }

    this->GetModelValue(curIndex,modelValue);
    if (isZero(modelValue))
    {
        for (unsigned int i = 0;i < 3;++i)
            newPoint[i] = oldPoint[i] + (k1[i] + 2.0 * k2[i] + 2.0 * k3[i]) / 5.0;

        return;
    }

    k4 = this->GetNextDirection(k3,modelValue,is2d,threadId);
    for (unsigned int i = 0;i < 3;++i)
        newPoint[i] = oldPoint[i] + (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) * m_StepProgression / 6.0;
}

bool BaseTractographyImageFilter::CheckIndexInImageBounds(IndexType &index, ImageBaseType *testImage)
{
    RegionType imageRegion = testImage->GetLargestPossibleRegion();
    for (unsigned int i = 0;i < ImageBaseType::ImageDimension;++i)
    {
        if ((index[i] < imageRegion.GetIndex()[i]) ||
            (index[i] >= imageRegion.GetIndex()[i] + imageRegion.GetSize()[i]))
            return false;
    }

    return true;
}

bool BaseTractographyImageFilter::isZero(VectorType &value)
{
    for (unsigned int i = 0;i < value.GetSize();++i)
    {
        if (value[i] != 0)
            return false;
    }

    return true;
}

} // end of namespace anima
