#include "animaBaseTractographyImageFilter.h"

#include <itkImageRegionIteratorWithIndex.h>
#include <itkMultiThreader.h>
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
    m_MinimalModelWeight = 0.25;

    m_ComputeLocalColors = true;
    m_HighestProcessedSeed = 0;
    m_ProgressReport = 0;
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

    unsigned int numSteps = std::floor(m_PointsToProcess.size() / (float)stepData);
    if (m_PointsToProcess.size() % stepData != 0)
        numSteps++;

    m_ProgressReport = new itk::ProgressReporter(this,0,numSteps);

    std::vector < FiberType > resultFibers;
    
    trackerArguments tmpStr;
    tmpStr.trackerPtr = this;
    for (unsigned int i = 0;i < this->GetNumberOfThreads();++i)
        tmpStr.resultFibersFromThreads.push_back(resultFibers);

    this->GetMultiThreader()->SetNumberOfThreads(this->GetNumberOfThreads());
    this->GetMultiThreader()->SetSingleMethod(this->ThreadTracker,&tmpStr);
    this->GetMultiThreader()->SingleMethodExecute();
    
    for (unsigned int j = 0;j < this->GetNumberOfThreads();++j)
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
    
    MaskImageIteratorType maskItr(m_MaskImage, m_InputImage->GetLargestPossibleRegion());
    m_PointsToProcess.clear();
    
    IndexType tmpIndex;
    PointType tmpPoint;
    ContinuousIndexType realIndex;
    
    m_FilteringValues.clear();
    double startN = -0.5 + 1.0 / (2.0 * m_NumberOfFibersPerPixel);
    double stepN = 1.0 / m_NumberOfFibersPerPixel;
    
    bool is2d = (m_InputImage->GetLargestPossibleRegion().GetSize()[2] <= 1);
    FiberType tmpFiber(1);
    
    while (!maskItr.IsAtEnd())
    {
        if (maskItr.Get() == 0)
        {
            ++maskItr;
            continue;
        }
        
        bool isAlreadyIn = false;
        for (unsigned int i = 0;i < m_FilteringValues.size();++i)
        {
            if (m_FilteringValues[i] == maskItr.Get())
            {
                isAlreadyIn = true;
                break;
            }
        }
        
        if (!isAlreadyIn)
            m_FilteringValues.push_back(maskItr.Get());
        
        tmpIndex = maskItr.GetIndex();
        
        if (is2d)
        {
            realIndex[2] = tmpIndex[2];
            for (unsigned int j = 0;j < m_NumberOfFibersPerPixel;++j)
            {
                realIndex[1] = tmpIndex[1] + startN + j * stepN;
                for (unsigned int i = 0;i < m_NumberOfFibersPerPixel;++i)
                {
                    realIndex[0] = tmpIndex[0] + startN + i * stepN;
                    m_MaskImage->TransformContinuousIndexToPhysicalPoint(realIndex,tmpPoint);
                    tmpFiber[0] = tmpPoint;
                    m_PointsToProcess.push_back(std::pair<FiberProgressType,FiberType>(Both,tmpFiber));
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
                        m_MaskImage->TransformContinuousIndexToPhysicalPoint(realIndex,tmpPoint);
                        tmpFiber[0] = tmpPoint;
                        m_PointsToProcess.push_back(std::pair<FiberProgressType,FiberType>(Both,tmpFiber));
                    }
                }
            }
        }
        
        ++maskItr;
    }
    
    std::cout << "Generated " << m_PointsToProcess.size() << " seed points from ROI mask" << std::endl;
}

ITK_THREAD_RETURN_TYPE BaseTractographyImageFilter::ThreadTracker(void *arg)
{
    itk::MultiThreader::ThreadInfoStruct *threadArgs = (itk::MultiThreader::ThreadInfoStruct *)arg;
    unsigned int nbThread = threadArgs->ThreadID;
    
    trackerArguments *tmpArg = (trackerArguments *)threadArgs->UserData;
    tmpArg->trackerPtr->ThreadTrack(nbThread,tmpArg->resultFibersFromThreads[nbThread]);
    
    return NULL;
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
        m_LockHighestProcessedSeed.Lock();

        if (m_HighestProcessedSeed >= highestToleratedSeedIndex)
        {
            m_LockHighestProcessedSeed.Unlock();
            continueLoop = false;
            continue;
        }

        unsigned int startPoint = m_HighestProcessedSeed;
        unsigned int endPoint = m_HighestProcessedSeed + stepData;
        if (endPoint > highestToleratedSeedIndex)
            endPoint = highestToleratedSeedIndex;

        m_HighestProcessedSeed = endPoint;

        m_LockHighestProcessedSeed.Unlock();

        this->ThreadedTrackComputer(numThread,resultFibers,startPoint,endPoint);

        m_LockHighestProcessedSeed.Lock();
        m_ProgressReport->CompletedPixel();
        m_LockHighestProcessedSeed.Unlock();
    }
}

void BaseTractographyImageFilter::ThreadedTrackComputer(unsigned int numThread, std::vector <FiberType> &resultFibers,
                                                        unsigned int startSeedIndex, unsigned int endSeedIndex)
{    
    FiberType tmpFiber;
    for (unsigned int i = startSeedIndex;i < endSeedIndex;++i)
    {
        tmpFiber = m_PointsToProcess[i].second;
        this->ComputeFiber(tmpFiber,m_PointsToProcess[i].first,numThread);
        
        if (tmpFiber.size() > m_MinLengthFiber / m_StepProgression)
            resultFibers.push_back(tmpFiber);
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
                m_MaskImage->TransformPhysicalPointToIndex(tmpPoint,tmpIndex);
                
                unsigned int maskValue = m_MaskImage->GetPixel(tmpIndex);
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
    vtkSmartPointer <vtkUnsignedCharArray> myColors = vtkUnsignedCharArray::New();
    vtkSmartPointer <vtkUnsignedCharArray> myCellColors = vtkUnsignedCharArray::New();
    myColors->SetNumberOfComponents (3);
    myCellColors->SetNumberOfComponents (3);
    PointType tmpDiff;
    
    for (unsigned int i = 0;i < filteredFibers.size();++i)
    {
        unsigned int npts = filteredFibers[i].size();
        vtkIdType* ids = new vtkIdType[npts];
        
        for (unsigned int j = 0;j < npts;++j)
        {
            ids[j] = myPoints->InsertNextPoint(filteredFibers[i][j][0],filteredFibers[i][j][1],filteredFibers[i][j][2]);

            if (m_ComputeLocalColors)
            {
                if ((j > 0)&&(j < npts-1))
                {
                    for(unsigned int k = 0;k < 3;++k)
                        tmpDiff[k] = filteredFibers[i][j+1][k] - filteredFibers[i][j-1][k];
                }
                else if (j == 0)
                {
                    for(unsigned int k = 0;k < 3;++k)
                        tmpDiff[k] = filteredFibers[i][j+1][k] - filteredFibers[i][j][k];
                }
                else
                {
                    for(unsigned int k = 0;k < 3;++k)
                        tmpDiff[k] = filteredFibers[i][j][k] - filteredFibers[i][j-1][k];
                }

                double normDiff = std::sqrt(tmpDiff[0] * tmpDiff[0] + tmpDiff[1] * tmpDiff[1] + tmpDiff[2] * tmpDiff[2]);

                for( unsigned int k = 0;k < 3;++k)
                {
                    double c = std::abs (tmpDiff[k] / normDiff) * 255.0;
                    myColors->InsertNextValue( (unsigned char)(c > 255.0 ? 255.0 : c) );
                }
            }
        }

        // cell color
        double normDiff = 0;
        for (unsigned int k = 0;k < 3;++k)
        {
            tmpDiff[k] = filteredFibers[i][0][k] - filteredFibers[i][npts-1][k];
            normDiff += tmpDiff[k] * tmpDiff[k];
        }

        normDiff = std::sqrt(normDiff);
        for (unsigned int i = 0;i < 3;++i)
        {
            double c = std::abs (tmpDiff[i] / normDiff) * 255.0;
            myCellColors->InsertNextValue ((unsigned char)(c > 255.0 ? 255.0 : c));
        }
        
        m_Output->InsertNextCell (VTK_POLY_LINE, npts, ids);
        delete[] ids;
    }
    
    m_Output->SetPoints (myPoints);
    if (m_ComputeLocalColors)
    {
        m_Output->GetPointData()->SetScalars (myColors);
        m_Output->GetCellData()->SetScalars (myCellColors);
        this->ComputeAdditionalScalarMaps();
    }
}

BaseTractographyImageFilter::FiberProcessVectorType
BaseTractographyImageFilter::ComputeFiber(BaseTractographyImageFilter::FiberType &fiber,
                                          BaseTractographyImageFilter::FiberProgressType ways,
                                          itk::ThreadIdType threadId)
{
    FiberProcessVectorType resVal;
    bool is2d = m_InputImage->GetLargestPossibleRegion().GetSize()[2] == 1;
    
    if ((ways == Forward)||(ways == Both))
    {
        PointType curPoint = fiber[fiber.size() - 1];
        PointType oldDir, newDir;
        
        ContinuousIndexType curIndex;
        IndexType curNearestIndex;
        bool continueLoop = true;
        m_MaskImage->TransformPhysicalPointToContinuousIndex(curPoint,curIndex);
        m_MaskImage->TransformPhysicalPointToIndex(curPoint,curNearestIndex);

        VectorType modelValue;
        modelValue.SetSize(1);
        
        while (continueLoop)
        {
            if (!m_CutMaskImage.IsNull())
            {
                if (m_CutMaskImage->GetPixel(curNearestIndex) != 0)
                {
                    continueLoop = false;
                    continue;
                }
            }
            
            if (!m_ForbiddenMaskImage.IsNull())
            {
                if (m_ForbiddenMaskImage->GetPixel(curNearestIndex) != 0)
                {
                    fiber.clear();
                    return resVal;
                }
            }
            
            this->GetModelValue(curIndex,modelValue);

            if (isZero(modelValue))
            {
                continueLoop = false;
                continue;
            }

            if (fiber.size() > m_MaxLengthFiber / m_StepProgression)
            {
                fiber.clear();
                return resVal;
            }
            
            if (fiber.size() <= 1)
                oldDir = this->GetModelPrincipalDirection(modelValue,is2d,threadId);
            else
            {
                for (unsigned int i = 0;i < 3;++i)
                    oldDir[i] = fiber[fiber.size() - 1][i] - fiber[fiber.size() - 2][i];
            }
            
            if (!this->CheckModelCompatibility(modelValue,threadId))
            {
                continueLoop = false;
                continue;
            }

            // Do some stuff to progress one step
            newDir = this->GetNextDirection(oldDir,modelValue,is2d,threadId);
            
            if (anima::ComputeOrientationAngle(oldDir, newDir) > m_MaxFiberAngle)
            {
                continueLoop = false;
                continue;
            }
            
            if (anima::ComputeScalarProduct(oldDir, newDir) < 0)
                anima::Revert(newDir,newDir);
            
            for (unsigned int i = 0;i < 3;++i)
                curPoint[i] = fiber[fiber.size() - 1][i] + m_StepProgression * newDir[i];
            
            m_MaskImage->TransformPhysicalPointToContinuousIndex(curPoint,curIndex);
            m_MaskImage->TransformPhysicalPointToIndex(curPoint,curNearestIndex);
            
            if (!this->CheckIndexInImageBounds(curIndex))
            {
                continueLoop = false;
                break;
            }
            
            // Add new point to fiber
            fiber.push_back(curPoint);
        }
    }

    if ((ways == Backward)||(ways == Both))
    {
        PointType curPoint = fiber[0];
        PointType oldDir, newDir;
        
        ContinuousIndexType curIndex;
        IndexType curNearestIndex;
        bool continueLoop = true;
        m_MaskImage->TransformPhysicalPointToContinuousIndex(curPoint,curIndex);
        m_MaskImage->TransformPhysicalPointToIndex(curPoint,curNearestIndex);
        
        VectorType modelValue;
        modelValue.SetSize(1);

        while (continueLoop)
        {
            if (!m_CutMaskImage.IsNull())
            {
                if (m_CutMaskImage->GetPixel(curNearestIndex) != 0)
                {
                    continueLoop = false;
                    continue;
                }
            }
            
            if (!m_ForbiddenMaskImage.IsNull())
            {
                if (m_ForbiddenMaskImage->GetPixel(curNearestIndex) != 0)
                {
                    fiber.clear();
                    return resVal;
                }
            }
            
            if (fiber.size() > m_MaxLengthFiber / m_StepProgression)
            {
                fiber.clear();
                return resVal;
            }
            
            this->GetModelValue(curIndex,modelValue);

            if (isZero(modelValue))
            {
                continueLoop = false;
                continue;
            }

            if (fiber.size() <= 1)
            {
                oldDir = this->GetModelPrincipalDirection(modelValue,is2d,threadId);
                for (unsigned int i = 0;i < 3;++i)
                    oldDir[i] *= -1.0;
            }
            else
            {
                for (unsigned int i = 0;i < 3;++i)
                    oldDir[i] = fiber[0][i] - fiber[1][i];
            }
            
            if (!this->CheckModelCompatibility(modelValue,threadId))
            {
                continueLoop = false;
                continue;
            }
            
            // Do some stuff to progress one step
            newDir = this->GetNextDirection(oldDir,modelValue,is2d,threadId);
            
            if (anima::ComputeOrientationAngle(oldDir, newDir) > m_MaxFiberAngle)
            {
                continueLoop = false;
                continue;
            }
            
            if (anima::ComputeScalarProduct(oldDir, newDir) < 0)
                anima::Revert(newDir,newDir);
            
            for (unsigned int i = 0;i < 3;++i)
                curPoint[i] = fiber[0][i] + m_StepProgression * newDir[i];
            
            m_MaskImage->TransformPhysicalPointToContinuousIndex(curPoint,curIndex);
            m_MaskImage->TransformPhysicalPointToIndex(curPoint,curNearestIndex);
            
            if (!this->CheckIndexInImageBounds(curIndex))
            {
                continueLoop = false;
                break;
            }
            
            // Now add point to fiber
            fiber.insert(fiber.begin(),1,curPoint);
        }
    }
    
    return resVal;
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
