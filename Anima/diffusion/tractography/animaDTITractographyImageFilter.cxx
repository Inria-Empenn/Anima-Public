#include "animaDTITractographyImageFilter.h"

#include <itkImageRegionIteratorWithIndex.h>
#include <itkSymmetricEigenAnalysis.h>
#include <itkMultiThreader.h>
#include <itkProgressReporter.h>

#include <animaVectorOperations.h>

#include <vtkPointData.h>

namespace anima
{

dtiTractographyImageFilter::dtiTractographyImageFilter()
{
    m_PointsToProcess.clear();
    
    m_NumberOfThreads = 1;
    m_NumberOfFibersPerPixel = 1;
    
    m_StepProgression = 1;
    m_StopFAThreshold = 0.1;
    m_MaxFiberAngle = M_PI / 2.0;

    m_ComputeLocalColors = true;
}

dtiTractographyImageFilter::~dtiTractographyImageFilter()
{
}

void dtiTractographyImageFilter::Update()
{
    this->PrepareTractography();
    m_Output = vtkPolyData::New();
    
    std::vector < FiberType > resultFibers;
    
    itk::MultiThreader::Pointer threadWorker = itk::MultiThreader::New();
    trackerArguments *tmpStr = new trackerArguments;
    tmpStr->trackerPtr = this;
    for (unsigned int i = 0;i < m_NumberOfThreads;++i)
        tmpStr->resultFibersFromThreads.push_back(resultFibers);
    
    threadWorker->SetNumberOfThreads(m_NumberOfThreads);
    threadWorker->SetSingleMethod(this->ThreadTracker,tmpStr);
    threadWorker->SingleMethodExecute();
    
    for (unsigned int j = 0;j < m_NumberOfThreads;++j)
    {
        resultFibers.insert(resultFibers.end(),tmpStr->resultFibersFromThreads[j].begin(),
                            tmpStr->resultFibersFromThreads[j].end());
    }
    
    delete tmpStr;
    
    std::cout << "\nTracked a total of " << resultFibers.size() << " fibers" << std::endl;
    resultFibers = this->FilterOutputFibers(resultFibers);
    std::cout << "Kept " << resultFibers.size() << " fibers after filtering" << std::endl;
    this->createVTKOutput(resultFibers);
}

void dtiTractographyImageFilter::PrepareTractography()
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

ITK_THREAD_RETURN_TYPE dtiTractographyImageFilter::ThreadTracker(void *arg)
{
    itk::MultiThreader::ThreadInfoStruct *threadArgs = (itk::MultiThreader::ThreadInfoStruct *)arg;
    unsigned int nbThread = threadArgs->ThreadID;
    
    trackerArguments *tmpArg = (trackerArguments *)threadArgs->UserData;
    tmpArg->trackerPtr->ThreadTrack(nbThread,tmpArg->resultFibersFromThreads[nbThread]);
    
    return NULL;
}

void dtiTractographyImageFilter::ThreadTrack(unsigned int numThread, std::vector <FiberType> &resultFibers)
{
    unsigned int dataSize = m_PointsToProcess.size();
    unsigned int numPoints = floor((double) dataSize / m_NumberOfThreads);
    
    unsigned int startPoint = numPoints * numThread;
    unsigned int endPoint = startPoint + numPoints;
    
    if (numThread+1 == m_NumberOfThreads)
        endPoint = dataSize;
    
    DTIInterpolatorPointer dtiInterpolator = DTIInterpolatorType::New();
    dtiInterpolator->SetInputImage(m_InputImage);
    
    // support progress methods/callbacks
    itk::ProgressReporter progress(this, numThread, endPoint - startPoint + 1);

    FiberType tmpFiber;
    for (unsigned int i = startPoint;i < endPoint;++i)
    {
        tmpFiber = m_PointsToProcess[i].second;
        this->ComputeFiber(tmpFiber,dtiInterpolator,m_PointsToProcess[i].first);
        
        if (tmpFiber.size() > m_MinLengthFiber / m_StepProgression)
            resultFibers.push_back(tmpFiber);

        progress.CompletedPixel();
    }
}

std::vector < std::vector <dtiTractographyImageFilter::PointType> >
dtiTractographyImageFilter::FilterOutputFibers(std::vector < FiberType > &fibers)
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

void dtiTractographyImageFilter::createVTKOutput(std::vector < FiberType > &filteredFibers)
{
    m_Output = vtkPolyData::New();
    m_Output->Initialize();
    m_Output->Allocate();
    
    vtkPoints* myPoints = vtkPoints::New();
    vtkUnsignedCharArray* myColors = vtkUnsignedCharArray::New();
    myColors->SetNumberOfComponents (3);
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

                double normDiff = sqrt(tmpDiff[0] * tmpDiff[0] + tmpDiff[1] * tmpDiff[1] + tmpDiff[2] * tmpDiff[2]);

                for( unsigned int k = 0;k < 3;++k)
                {
                    double c = fabs (tmpDiff[k] / normDiff) * 255.0;
                    myColors->InsertNextValue( (unsigned char)(c > 255.0 ? 255.0 : c) );
                }
            }
        }
        
        m_Output->InsertNextCell (VTK_POLY_LINE, npts, ids);
        delete[] ids;
    }
    
    m_Output->SetPoints (myPoints);
    if (m_ComputeLocalColors)
        m_Output->GetPointData()->SetScalars (myColors);
    myPoints->Delete();
    myColors->Delete();
}

dtiTractographyImageFilter::FiberProcessVectorType
dtiTractographyImageFilter::ComputeFiber(dtiTractographyImageFilter::FiberType &fiber,
                                         dtiTractographyImageFilter::DTIInterpolatorType * dtiInterpolator,
                                         dtiTractographyImageFilter::FiberProgressType ways)
{
    FiberProcessVectorType resVal;
    bool is2d = m_InputImage->GetLargestPossibleRegion().GetSize()[2] == 1;
    
    if ((ways == Forward)||(ways == Both))
    {
        PointType curPoint = fiber[fiber.size() - 1];
        std::vector <double> oldDir(3,0), newDir(3,0);
        
        ContinuousIndexType curIndex;
        IndexType curNearestIndex;
        bool continueLoop = true;
        m_MaskImage->TransformPhysicalPointToContinuousIndex(curPoint,curIndex);
        m_MaskImage->TransformPhysicalPointToIndex(curPoint,curNearestIndex);
        
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
            
            VectorType tmpDTIValue = dtiInterpolator->EvaluateAtContinuousIndex(curIndex);
            
            if (fiber.size() > m_MaxLengthFiber / m_StepProgression)
            {
                fiber.clear();
                return resVal;
            }
            
            if (fiber.size() <= 1)
                oldDir = this->GetDTIPrincipalDirection(tmpDTIValue,is2d);
            else
            {
                for (unsigned int i = 0;i < 3;++i)
                    oldDir[i] = fiber[fiber.size() - 1][i] - fiber[fiber.size() - 2][i];
            }
            
            double fa = this->GetFA(tmpDTIValue);
            
            if (fa < m_StopFAThreshold)
            {
                continueLoop = false;
                continue;
            }
            
            // Do some stuff to progress one step
            newDir = this->GetDTIPrincipalDirection(tmpDTIValue,is2d);
            
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
            
            if (!dtiInterpolator->IsInsideBuffer(curIndex))
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
        std::vector <double> oldDir(3,0), newDir(3,0);
        
        ContinuousIndexType curIndex;
        IndexType curNearestIndex;
        bool continueLoop = true;
        m_MaskImage->TransformPhysicalPointToContinuousIndex(curPoint,curIndex);
        m_MaskImage->TransformPhysicalPointToIndex(curPoint,curNearestIndex);
        
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
            
            VectorType tmpDTIValue = dtiInterpolator->EvaluateAtContinuousIndex(curIndex);
            
            if (fiber.size() <= 1)
            {
                oldDir = this->GetDTIPrincipalDirection(tmpDTIValue,is2d);
                for (unsigned int i = 0;i < 3;++i)
                    oldDir[i] *= -1.0;
            }
            else
            {
                for (unsigned int i = 0;i < 3;++i)
                    oldDir[i] = fiber[0][i] - fiber[1][i];
            }
            
            double fa = this->GetFA(tmpDTIValue);
            
            if (fa < m_StopFAThreshold)
            {
                continueLoop = false;
                continue;
            }
            
            // Do some stuff to progress one step
            newDir = this->GetDTIPrincipalDirection(tmpDTIValue,is2d);
            
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
            
            if (!dtiInterpolator->IsInsideBuffer(curIndex))
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

double dtiTractographyImageFilter::GetFA(VectorType &dtiLogValue)
{
    typedef vnl_matrix <double> MatrixType;
    
    itk::SymmetricEigenAnalysis <MatrixType,vnl_diag_matrix <double>,MatrixType> EigenAnalysis(3);
    
    vnl_matrix <double> tmpMat(3,3);
    
    unsigned int pos = 0;
    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j <= i;++j)
        {
            tmpMat(i,j) = dtiLogValue[pos];
            if (j != i)
                tmpMat(j,i) = tmpMat(i,j);
            ++pos;
        }
    
    vnl_diag_matrix <double> eVals(3);
    EigenAnalysis.ComputeEigenValues(tmpMat,eVals);
    
    double meanEvals = 0;
    for (unsigned int i = 0;i < 3;++i)
    {
        eVals[i] = exp(eVals[i]);
        meanEvals += eVals[i];
    }
    
    meanEvals /= 3.0;
    
    double num = 0;
    double denom = 0;
    for (unsigned int i = 0;i < 3;++i)
    {
        num += (eVals[i] - meanEvals) * (eVals[i] - meanEvals);
        denom += eVals[i] * eVals[i];
    }
    
    return sqrt(num / (2 * denom));
}

std::vector <double> dtiTractographyImageFilter::GetDTIPrincipalDirection(VectorType &dtiLogValue, bool is2d)
{
    typedef vnl_matrix <double> MatrixType;
    
    itk::SymmetricEigenAnalysis <MatrixType,vnl_diag_matrix <double>,MatrixType> EigenAnalysis(3);
    
    MatrixType tmpMat(3,3);
    vnl_diag_matrix <double> eVals(3);
    MatrixType eVec(3,3);
    
    unsigned int pos = 0;
    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j <= i;++j)
        {
            tmpMat(i,j) = dtiLogValue[pos];
            if (j != i)
                tmpMat(j,i) = tmpMat(i,j);
            ++pos;
        }
    
    EigenAnalysis.ComputeEigenValuesAndVectors(tmpMat,eVals,eVec);
    std::vector <double> resDir(3,0);
    
    for (unsigned int i = 0;i < 3;++i)
        resDir[i] = eVec(2,i);
    
    if (is2d)
    {
        resDir[2] = 0;
        
        double norm = 0;
        for (unsigned int i = 0;i < 2;++i)
            norm += resDir[i] * resDir[i];
        norm = sqrt(norm);
        
        for (unsigned int i = 0;i < 2;++i)
            resDir[i] /= norm;
    }
    
    return resDir;
}

} // end of namespace anima
