#include "animaBaseProbabilisticTractographyImageFilter.h"

#include <itkImageRegionIteratorWithIndex.h>
#include <itkLinearInterpolateImageFunction.h>

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkExtractImageFilter.h>
#include <itkImageMomentsCalculator.h>

#include <animaVectorOperations.h>
#include <animaLogarithmFunctions.h>

#include <vnl/algo/vnl_matrix_inverse.h>

#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>

#include <animaKMeansFilter.h>

#include <ctime>

namespace anima
{

template <class TInputModelImageType>
BaseProbabilisticTractographyImageFilter <TInputModelImageType>
::BaseProbabilisticTractographyImageFilter()
{
    m_PointsToProcess.clear();

    m_NumberOfFibersPerPixel = 1;
    m_NumberOfParticles = 1000;
    m_MinimalNumberOfParticlesPerClass = 10;

    m_ResamplingThreshold = 0.8;

    m_StepProgression = 1.0;

    m_KappaOfPriorDistribution = 30.0;

    m_MinLengthFiber = 10.0;
    m_MaxLengthFiber = 150.0;

    m_PositionDistanceFuseThreshold = 0.5;
    m_KappaSplitThreshold = 30.0;

    m_ClusterDistance = 1;

    m_ComputeLocalColors = true;
    m_MAPMergeFibers = true;

    m_InitialColinearityDirection = Center;
    m_InitialDirectionMode = Weight;

    m_Generators.clear();

    m_HighestProcessedSeed = 0;
    m_ProgressReport = 0;
}

template <class TInputModelImageType>
BaseProbabilisticTractographyImageFilter <TInputModelImageType>
::~BaseProbabilisticTractographyImageFilter()
{
    if (m_ProgressReport)
        delete m_ProgressReport;
}

template <class TInputModelImageType>
void
BaseProbabilisticTractographyImageFilter <TInputModelImageType>
::Update()
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

    FiberProcessVectorType resultFibers;
    ListType resultWeights;

    trackerArguments tmpStr;
    tmpStr.trackerPtr = this;
    tmpStr.resultFibersFromThreads.resize(this->GetNumberOfWorkUnits());
    tmpStr.resultWeightsFromThreads.resize(this->GetNumberOfWorkUnits());

    for (unsigned int i = 0;i < this->GetNumberOfWorkUnits();++i)
    {
        tmpStr.resultFibersFromThreads[i] = resultFibers;
        tmpStr.resultWeightsFromThreads[i] = resultWeights;
    }

    this->GetMultiThreader()->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    this->GetMultiThreader()->SetSingleMethod(this->ThreadTracker,&tmpStr);
    this->GetMultiThreader()->SingleMethodExecute();

    for (unsigned int j = 0;j < this->GetNumberOfWorkUnits();++j)
    {
        resultFibers.insert(resultFibers.end(),tmpStr.resultFibersFromThreads[j].begin(),tmpStr.resultFibersFromThreads[j].end());
        resultWeights.insert(resultWeights.end(),tmpStr.resultWeightsFromThreads[j].begin(),tmpStr.resultWeightsFromThreads[j].end());
    }

    std::cout << "\nKept " << resultFibers.size() << " fibers after filtering" << std::endl;
    this->createVTKOutput(resultFibers, resultWeights);
}

template <class TInputModelImageType>
void
BaseProbabilisticTractographyImageFilter <TInputModelImageType>
::PrepareTractography()
{
    if (!m_B0Image)
        itkExceptionMacro("No B0 image, required");

    if (!m_NoiseImage)
        itkExceptionMacro("No sigma square noise image, required");

    m_B0Interpolator = ScalarInterpolatorType::New();
    m_B0Interpolator->SetInputImage(m_B0Image);

    m_NoiseInterpolator = ScalarInterpolatorType::New();
    m_NoiseInterpolator->SetInputImage(m_NoiseImage);

    // Initialize random generator
    m_Generators.resize(this->GetNumberOfWorkUnits());

    std::mt19937 motherGenerator(time(0));

    for (unsigned int i = 0;i < this->GetNumberOfWorkUnits();++i)
        m_Generators[i] = std::mt19937(motherGenerator());

    bool is2d = m_InputModelImage->GetLargestPossibleRegion().GetSize()[2] == 1;
    if (is2d && (m_InitialColinearityDirection == Top))
        m_InitialColinearityDirection = Front;
    if (is2d && (m_InitialColinearityDirection == Bottom))
        m_InitialColinearityDirection = Back;

    // If needed, ensure DWI gravity center is computed
    if ((m_InitialColinearityDirection == Outward)||(m_InitialColinearityDirection == Center))
    {
        itk::ImageMomentsCalculator <ScalarImageType>::Pointer momentsCalculator = itk::ImageMomentsCalculator <ScalarImageType>::New();
        momentsCalculator->SetImage(m_B0Image);
        momentsCalculator->Compute();
        m_DWIGravityCenter = momentsCalculator->GetCenterOfGravity();
    }

    typedef itk::ImageRegionIteratorWithIndex <MaskImageType> MaskImageIteratorType;

    MaskImageIteratorType maskItr(m_SeedMask, m_InputModelImage->GetLargestPossibleRegion());
    m_PointsToProcess.clear();

    IndexType tmpIndex;
    PointType tmpPoint;
    ContinuousIndexType realIndex;

    m_FilteringValues.clear();
    double startN = -0.5 + 1.0 / (2.0 * m_NumberOfFibersPerPixel);
    double stepN = 1.0 / m_NumberOfFibersPerPixel;
    FiberType tmpFiber(1);

    if (m_FilterMask)
    {
        MaskImageIteratorType filterItr(m_FilterMask, m_InputModelImage->GetLargestPossibleRegion());
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

    while (!maskItr.IsAtEnd())
    {
        if (maskItr.Get() == 0)
        {
            ++maskItr;
            continue;
        }

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
                    m_SeedMask->TransformContinuousIndexToPhysicalPoint(realIndex,tmpPoint);
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

                        m_SeedMask->TransformContinuousIndexToPhysicalPoint(realIndex,tmpPoint);
                        tmpFiber[0] = tmpPoint;
                        m_PointsToProcess.push_back(tmpFiber);
                    }
                }
            }
        }

        ++maskItr;
    }

    std::cout << "Generated " << m_PointsToProcess.size() << " seed points from ROI mask" << std::endl;
}

template <class TInputModelImageType>
typename BaseProbabilisticTractographyImageFilter <TInputModelImageType>::InterpolatorType *
BaseProbabilisticTractographyImageFilter <TInputModelImageType>::GetModelInterpolator()
{
    typedef itk::LinearInterpolateImageFunction <InputModelImageType> InternalInterpolatorType;

    typename InternalInterpolatorType::Pointer outInterpolator = InternalInterpolatorType::New();
    outInterpolator->SetInputImage(m_InputModelImage);

    outInterpolator->Register();
    return outInterpolator;
}

template <class TInputModelImageType>
ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION
BaseProbabilisticTractographyImageFilter <TInputModelImageType>
::ThreadTracker(void *arg)
{
    itk::MultiThreaderBase::WorkUnitInfo *threadArgs = (itk::MultiThreaderBase::WorkUnitInfo *)arg;
    unsigned int nbThread = threadArgs->WorkUnitID;

    trackerArguments *tmpArg = (trackerArguments *)threadArgs->UserData;
    tmpArg->trackerPtr->ThreadTrack(nbThread,tmpArg->resultFibersFromThreads[nbThread],tmpArg->resultWeightsFromThreads[nbThread]);

    return ITK_THREAD_RETURN_DEFAULT_VALUE;
}

template <class TInputModelImageType>
void
BaseProbabilisticTractographyImageFilter <TInputModelImageType>
::ThreadTrack(unsigned int numThread, FiberProcessVectorType &resultFibers,
              ListType &resultWeights)
{
    bool continueLoop = true;
    unsigned int highestToleratedSeedIndex = m_PointsToProcess.size();

    unsigned int stepData = std::min((int)m_PointsToProcess.size(),100);
    if (stepData == 0)
        stepData = 1;

    while (continueLoop)
    {
        m_LockHighestProcessedSeed.lock();

        if (m_HighestProcessedSeed >= highestToleratedSeedIndex)
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

        this->ThreadedTrackComputer(numThread,resultFibers,resultWeights,startPoint,endPoint);

        m_LockHighestProcessedSeed.lock();
        m_ProgressReport->CompletedPixel();
        m_LockHighestProcessedSeed.unlock();
    }
}

template <class TInputModelImageType>
void
BaseProbabilisticTractographyImageFilter <TInputModelImageType>
::ThreadedTrackComputer(unsigned int numThread, FiberProcessVectorType &resultFibers,
                        ListType &resultWeights, unsigned int startSeedIndex,
                        unsigned int endSeedIndex)
{
    InterpolatorPointer modelInterpolator = this->GetModelInterpolator();
    FiberProcessVectorType tmpFibers;
    ListType tmpWeights;
    ContinuousIndexType startIndex;

    for (unsigned int i = startSeedIndex;i < endSeedIndex;++i)
    {
        m_SeedMask->TransformPhysicalPointToContinuousIndex(m_PointsToProcess[i][0],startIndex);

        tmpFibers = this->ComputeFiber(m_PointsToProcess[i], modelInterpolator, numThread, tmpWeights);

        tmpFibers = this->FilterOutputFibers(tmpFibers, tmpWeights);

        for (unsigned int j = 0;j < tmpFibers.size();++j)
        {
            if (tmpFibers[j].size() > m_MinLengthFiber / m_StepProgression)
            {
                resultFibers.push_back(tmpFibers[j]);
                resultWeights.push_back(tmpWeights[j]);
            }
        }
    }
}

template <class TInputModelImageType>
typename BaseProbabilisticTractographyImageFilter <TInputModelImageType>::FiberProcessVectorType
BaseProbabilisticTractographyImageFilter <TInputModelImageType>
::FilterOutputFibers(FiberProcessVectorType &fibers, ListType &weights)
{
    FiberProcessVectorType resVal;
    ListType tmpWeights = weights;
    weights.clear();

    if ((m_FilteringValues.size() > 0)||(m_ForbiddenMask))
    {
        MembershipType touchingLabels;
        IndexType tmpIndex;
        PointType tmpPoint;

        for (unsigned int i = 0;i < fibers.size();++i)
        {
            touchingLabels.clear();
            bool forbiddenTouched = false;

            for (unsigned int j = 0;j < fibers[i].size();++j)
            {
                tmpPoint = fibers[i][j];
                m_SeedMask->TransformPhysicalPointToIndex(tmpPoint,tmpIndex);

                unsigned int maskValue = 0;
                unsigned int forbiddenMaskValue = 0;

                if (m_FilterMask)
                    maskValue = m_FilterMask->GetPixel(tmpIndex);

                if (m_ForbiddenMask)
                    forbiddenMaskValue = m_ForbiddenMask->GetPixel(tmpIndex);

                if (forbiddenMaskValue != 0)
                {
                    forbiddenTouched = true;
                    break;
                }

                if (maskValue != 0)
                {
                    bool alreadyIn = false;
                    for (unsigned int k = 0;k < touchingLabels.size();++k)
                    {
                        if (maskValue == touchingLabels[k])
                        {
                            alreadyIn = true;
                            break;
                        }
                    }

                    if (!alreadyIn)
                        touchingLabels.push_back(maskValue);
                }
            }

            if (forbiddenTouched)
                continue;

            if (touchingLabels.size() == m_FilteringValues.size())
            {
                resVal.push_back(fibers[i]);
                weights.push_back(tmpWeights[i]);
            }
        }
    }
    else
    {
        resVal = fibers;
        weights = tmpWeights;
    }

    return resVal;
}

template <class TInputModelImageType>
void
BaseProbabilisticTractographyImageFilter <TInputModelImageType>
::createVTKOutput(FiberProcessVectorType &filteredFibers, ListType &filteredWeights)
{
    m_Output = vtkPolyData::New();
    m_Output->Initialize();
    m_Output->Allocate();

    vtkSmartPointer <vtkPoints> myPoints = vtkPoints::New();
    vtkSmartPointer <vtkUnsignedCharArray> myColors = vtkUnsignedCharArray::New();
    vtkSmartPointer <vtkUnsignedCharArray> myCellColors = vtkUnsignedCharArray::New();
    myColors->SetNumberOfComponents (3);
    myCellColors->SetNumberOfComponents (3);

    vtkSmartPointer <vtkDoubleArray> weights = vtkDoubleArray::New();
    weights->SetNumberOfComponents(1);
    weights->SetName("Fiber weights");

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
                    double c = std::abs (tmpDiff[k] / normDiff) * 255.0;
                    myColors->InsertNextValue( (unsigned char)(c > 255.0 ? 255.0 : c) );
                }
            }

            weights->InsertNextValue( filteredWeights[i] );
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

    m_Output->SetPoints(myPoints);
    if (m_ComputeLocalColors)
    {
        m_Output->GetPointData()->SetScalars(myColors);
        m_Output->GetCellData()->SetScalars(myCellColors);
        this->ComputeAdditionalScalarMaps();
    }

    // Add particle weights to data
    m_Output->GetPointData()->AddArray(weights);
}

template <class TInputModelImageType>
typename BaseProbabilisticTractographyImageFilter <TInputModelImageType>::FiberProcessVectorType
BaseProbabilisticTractographyImageFilter <TInputModelImageType>
::ComputeFiber(FiberType &fiber, InterpolatorPointer &modelInterpolator,
               unsigned int numThread, ListType &resultWeights)
{
    unsigned int numberOfClasses = 1;

    FiberWorkType fiberComputationData;
    fiberComputationData.fiberParticles.resize(m_NumberOfParticles);
    for (unsigned int i = 0;i < m_NumberOfParticles;++i)
        fiberComputationData.fiberParticles[i] = fiber;

    fiberComputationData.particleWeights = ListType(m_NumberOfParticles, 1.0 / m_NumberOfParticles);
    fiberComputationData.stoppedParticles = std::vector <bool> (m_NumberOfParticles,false);
    fiberComputationData.classSizes = MembershipType(numberOfClasses,m_NumberOfParticles);
    fiberComputationData.classWeights = ListType(numberOfClasses,1.0 / numberOfClasses);
    //We need membership vectors in each direction
    fiberComputationData.classMemberships = MembershipType(m_NumberOfParticles,0);
    fiberComputationData.reverseClassMemberships.resize(numberOfClasses);
    MembershipType tmpVec(m_NumberOfParticles);
    for (unsigned int j = 0;j < m_NumberOfParticles;++j)
        tmpVec[j] = j;
    fiberComputationData.reverseClassMemberships[0] = tmpVec;

    ListType logWeightVals(m_NumberOfParticles, 1.0 / m_NumberOfParticles);
    ListType oldFiberWeights(m_NumberOfParticles, 1.0 / m_NumberOfParticles);
    std::vector <bool> emptyClasses;
    ListType tmpVector(m_NumberOfParticles,0);

    ListType logWeightSums(numberOfClasses,0);
    ListType effectiveNumberOfParticles(numberOfClasses,0);

    DirectionVectorType previousDirections(m_NumberOfParticles);

    // Data structures for resampling
    FiberProcessVectorType fiberParticlesCopy;
    DirectionVectorType previousDirectionsCopy;
    ListType weightSpecificClassValues;
    FiberProcessVectorType fiberTrash;
    std::vector <bool> usedFibers;

    // Here to constrain directions to 2D plane if needed
    bool is2d = m_InputModelImage->GetLargestPossibleRegion().GetSize()[2] == 1;

    VectorType modelValue(m_ModelDimension);

    Vector3DType sampling_direction(0.0), newDirection;
    PointType currentPoint;
    ContinuousIndexType currentIndex, newIndex;
    IndexType closestIndex;

    unsigned int numIter = 0;
    bool stopLoop = false;
    while (!stopLoop)
    {
        ++numIter;

        // Store previous weights for the resampling step
        if (numIter > 1)
            oldFiberWeights = fiberComputationData.particleWeights;

        logWeightSums.resize(numberOfClasses);
        std::fill(logWeightSums.begin(),logWeightSums.end(),0.0);
        for (unsigned int i = 0;i < m_NumberOfParticles;++i)
        {
            // Do not compute trashed fibers
            if (fiberComputationData.stoppedParticles[i])
                continue;

            currentPoint = fiberComputationData.fiberParticles[i].back();

            m_SeedMask->TransformPhysicalPointToContinuousIndex(currentPoint,currentIndex);

            // Trash fiber if it goes outside of the brain
            if (!modelInterpolator->IsInsideBuffer(currentIndex))
            {
                fiberComputationData.stoppedParticles[i] = true;
                fiberComputationData.particleWeights[i] = 0;
                continue;
            }

            // Trash fiber if it goes through the cut mask
            m_SeedMask->TransformPhysicalPointToIndex(currentPoint,closestIndex);

            if (m_CutMask)
            {
                if (m_CutMask->GetPixel(closestIndex) != 0)
                {
                    fiberComputationData.stoppedParticles[i] = true;
                    fiberComputationData.particleWeights[i] = 0;
                    continue;
                }
            }

            // Computes diffusion information at current position
            modelValue.Fill(0.0);
            double estimatedNoiseValue = 20.0;
            this->ComputeModelValue(modelInterpolator, currentIndex, modelValue);
            double estimatedB0Value = m_B0Interpolator->EvaluateAtContinuousIndex(currentIndex);
            estimatedNoiseValue = m_NoiseInterpolator->EvaluateAtContinuousIndex(currentIndex);

            if (!this->CheckModelProperties(estimatedB0Value,estimatedNoiseValue,modelValue,numThread))
            {
                fiberComputationData.stoppedParticles[i] = true;
                fiberComputationData.particleWeights[i] = 0;
                continue;
            }

            // Set initial direction to the principal eigenvector of the tensor
            if (numIter == 1)
            {
                Vector3DType initDir(0.0);
                switch (m_InitialColinearityDirection)
                {
                    case Top:
                        initDir[2] = 1;
                        break;
                    case Bottom:
                        initDir[2] = -1;
                        break;
                    case Left:
                        initDir[0] = -1;
                        break;
                    case Right:
                        initDir[0] = 1;
                        break;
                    case Front:
                        initDir[1] = -1;
                        break;
                    case Back:
                        initDir[1] = 1;
                        break;
                    case Outward:
                        for (unsigned int j = 0;j < InputModelImageType::ImageDimension;++j)
                            initDir[j] = currentPoint[j] - m_DWIGravityCenter[j];
                        break;
                    case Center:
                    default:
                        for (unsigned int j = 0;j < InputModelImageType::ImageDimension;++j)
                            initDir[j] = m_DWIGravityCenter[j] - currentPoint[j];
                        break;
                }

                if (is2d)
                    initDir[2] = 0;
                initDir.Normalize();

                previousDirections[i] = this->InitializeFirstIterationFromModel(initDir,modelValue,numThread);
            }

            // Propose a new direction based on the previous one and the diffusion information at current position
            double log_prior = 0, log_proposal = 0;
            newDirection = this->ProposeNewDirection(previousDirections[i], modelValue, sampling_direction, log_prior,
                                                     log_proposal, m_Generators[numThread], numThread);

            // Update the position of the particle
            for (unsigned int j = 0;j < InputModelImageType::ImageDimension;++j)
                currentPoint[j] += m_StepProgression * newDirection[j];

            // Log-weight update must be done at new position (except for prior and proposal)
            m_SeedMask->TransformPhysicalPointToContinuousIndex(currentPoint,newIndex);

            // Set the new proposed direction as the current direction
            previousDirections[i] = newDirection;

            modelValue.Fill(0.0);

            if (!modelInterpolator->IsInsideBuffer(newIndex))
            {
                fiberComputationData.stoppedParticles[i] = true;
                fiberComputationData.particleWeights[i] = 0;
                continue;
            }

            fiberComputationData.fiberParticles[i].push_back(currentPoint);

            this->ComputeModelValue(modelInterpolator, newIndex, modelValue);
            estimatedB0Value = m_B0Interpolator->EvaluateAtContinuousIndex(newIndex);
            estimatedNoiseValue = m_NoiseInterpolator->EvaluateAtContinuousIndex(newIndex);

            // Update the weight of the particle
            double updateWeightLogVal = this->ComputeLogWeightUpdate(estimatedB0Value, estimatedNoiseValue, newDirection,
                                                                     modelValue, log_prior, log_proposal, numThread);

            logWeightVals[i] = updateWeightLogVal + anima::safe_log(oldFiberWeights[i]);
        }

        // Continue only if some particles are still moving
        stopLoop = true;
        for (unsigned int i = 0;i < m_NumberOfParticles;++i)
        {
            if (!fiberComputationData.stoppedParticles[i])
            {
                stopLoop = false;
                break;
            }
        }

        if (stopLoop)
            continue;

        emptyClasses.resize(numberOfClasses);
        // Computes weight sum for further weight normalization
        for (unsigned int i = 0;i < numberOfClasses;++i)
        {
            emptyClasses[i] = false;
            unsigned int classSize = fiberComputationData.reverseClassMemberships[i].size();
            tmpVector.clear();

            for (unsigned int j = 0;j < classSize;++j)
            {
                if (!fiberComputationData.stoppedParticles[fiberComputationData.reverseClassMemberships[i][j]])
                    tmpVector.push_back(logWeightVals[fiberComputationData.reverseClassMemberships[i][j]]);
            }

            if (tmpVector.size() != 0)
                logWeightSums[i] = anima::ExponentialSum(tmpVector);
            else
            {
                logWeightSums[i] = 0;
                emptyClasses[i] = true;
            }
        }

        // Weight normalization
        tmpVector.clear();
        for (unsigned int i = 0;i < numberOfClasses;++i)
        {
            if (!emptyClasses[i])
                tmpVector.push_back(anima::safe_log(fiberComputationData.classWeights[i]) + logWeightSums[i]);
        }

        double tmpSum = 0;
        double logSumTmpVector = anima::ExponentialSum(tmpVector);

        for (unsigned int i = 0;i < numberOfClasses;++i)
        {
            if (!emptyClasses[i])
            {
                double t = std::exp(anima::safe_log(fiberComputationData.classWeights[i]) + logWeightSums[i] - logSumTmpVector);
                fiberComputationData.classWeights[i] = t;
                tmpSum += t;
            }
            else
                fiberComputationData.classWeights[i] = 0;
        }

        for (unsigned int i = 0;i < numberOfClasses;++i)
            fiberComputationData.classWeights[i] /= tmpSum;

        for (unsigned int i = 0;i < m_NumberOfParticles;++i)
        {
            if (fiberComputationData.stoppedParticles[i])
            {
                fiberComputationData.particleWeights[i] = 0;
                continue;
            }

            double tmpWeight = logWeightSums[fiberComputationData.classMemberships[i]];

            if (std::isfinite(tmpWeight))
                logWeightVals[i] -= tmpWeight;

            fiberComputationData.particleWeights[i] = std::exp(logWeightVals[i]);
            // Q: shouldn't we be treating this case as an empty cluster that shouldn't even exist?
        }

        // Resampling if necessary, done class by class
        effectiveNumberOfParticles.resize(numberOfClasses);
        std::fill(effectiveNumberOfParticles.begin(),effectiveNumberOfParticles.end(),0);
        for (unsigned int i = 0;i < m_NumberOfParticles;++i)
        {
            double weight = fiberComputationData.particleWeights[i];
            effectiveNumberOfParticles[fiberComputationData.classMemberships[i]] += weight * weight;
        }

        for (unsigned int m = 0;m < numberOfClasses;++m)
        {
            if (effectiveNumberOfParticles[m] != 0)
                effectiveNumberOfParticles[m] = 1.0 / effectiveNumberOfParticles[m];
            else
                continue; // Q: shouldn't we be treating this case as an empty cluster that shouldn't even exist? (same as previous Q)

            // Actual class resampling
            if (effectiveNumberOfParticles[m] < m_ResamplingThreshold * fiberComputationData.classSizes[m])
            {
                weightSpecificClassValues.resize(fiberComputationData.classSizes[m]);
                previousDirectionsCopy.resize(fiberComputationData.classSizes[m]);
                fiberParticlesCopy.resize(fiberComputationData.classSizes[m]);

                for (unsigned int i = 0;i < fiberComputationData.classSizes[m];++i)
                {
                    unsigned int posIndex = fiberComputationData.reverseClassMemberships[m][i];
                    weightSpecificClassValues[i] = fiberComputationData.particleWeights[posIndex];
                    previousDirectionsCopy[i] = previousDirections[posIndex];
                    fiberParticlesCopy[i] = fiberComputationData.fiberParticles[posIndex];
                }

                std::discrete_distribution<> dist(weightSpecificClassValues.begin(),weightSpecificClassValues.end());
                usedFibers.resize(fiberComputationData.classSizes[m]);
                std::fill(usedFibers.begin(),usedFibers.end(),false);

                for (unsigned int i = 0;i < fiberComputationData.classSizes[m];++i)
                {
                    unsigned int z = dist(m_Generators[numThread]);
                    unsigned int iReal = fiberComputationData.reverseClassMemberships[m][i];
                    previousDirections[iReal] = previousDirectionsCopy[z];
                    fiberComputationData.fiberParticles[iReal] = fiberParticlesCopy[z];
                    // In all of this, we suppose that stopped particles have zero weights and will therefore
                    // be lost when resampling
                    fiberComputationData.stoppedParticles[iReal] = false;
                    usedFibers[z] = true;
                }

                for (unsigned int i = 0;i < fiberComputationData.classSizes[m];++i)
                {
                    unsigned int iReal = fiberComputationData.reverseClassMemberships[m][i];

                    if (!usedFibers[i])
                    {
                        if (oldFiberWeights[iReal] > m_FiberTrashThreshold / m_NumberOfParticles)
                        {
                            if (fiberComputationData.particleWeights[iReal] != 0)
                                fiberParticlesCopy[i].pop_back();

                            // The fiber trash used to contain fibers that were lost with a sufficient weight
                            // However, using way too much memory so removed for now
                            //if (fiberParticlesCopy[i].size() > m_MinLengthFiber / m_StepProgression)
                            //    fiberTrash.push_back(fiberParticlesCopy[i]);
                        }
                    }
                }

                // Update only weightVals, oldWeightVals will get updated when starting back the loop
                // Same here for stopped fibers, they get rejected when resampling
                for (unsigned int i = 0;i < fiberComputationData.classSizes[m];++i)
                    fiberComputationData.particleWeights[fiberComputationData.reverseClassMemberships[m][i]] = 1.0 / fiberComputationData.classSizes[m];
            }
        }

        // We need stopping criterions
        // Length is easy, given that each step is constant we just need to check the fiber size: numIter
        // Example :
        if (numIter > m_MaxLengthFiber / m_StepProgression)
            stopLoop = true;

        numberOfClasses = this->UpdateClassesMemberships(fiberComputationData,previousDirections,m_Generators[numThread]);

        for (unsigned int i = 0;i < fiberComputationData.particleWeights.size();++i)
        {
            if (!std::isfinite(fiberComputationData.particleWeights[i]))
                itkExceptionMacro("Nan weights after update class membership");
        }
    }

    // Now that we're done, if we don't keep individual particles, merge them cluster by cluster
    if (m_MAPMergeFibers)
    {
        FiberProcessVectorType mergedOutput,classMergedOutput;
        for (unsigned int i = 0;i < numberOfClasses;++i)
        {
            this->MergeParticleClassFibers(fiberComputationData,classMergedOutput,i);
            mergedOutput.insert(mergedOutput.end(),classMergedOutput.begin(),classMergedOutput.end());
        }

        fiberComputationData.fiberParticles = mergedOutput;
    }

    if (m_MAPMergeFibers)
        resultWeights = fiberComputationData.classWeights;
    else
        resultWeights = fiberComputationData.particleWeights;

    return fiberComputationData.fiberParticles;
}

template <class TInputModelImageType>
unsigned int
BaseProbabilisticTractographyImageFilter <TInputModelImageType>
::UpdateClassesMemberships(FiberWorkType &fiberData, DirectionVectorType &directions, std::mt19937 &random_generator)
{
    const unsigned int p = PointType::PointDimension;
    typedef anima::KMeansFilter <PointType,p> KMeansFilterType;
    unsigned int numClasses = fiberData.classSizes.size();

    // Deciding on cluster merges
    FiberProcessVectorType mapMergedFibersRef, mapMergedFibersFlo;
    unsigned int newNumClasses = numClasses;

    MembershipType classesFusion(numClasses);
    for (unsigned int i = 0;i < numClasses;++i)
        classesFusion[i] = i;

    // As described in IPMI, we take the input classes and first try to fuse them
    // This is based on a range of possible criterions specified by the user
    if (numClasses > 1)
    {
        for (unsigned int i = 0;i < numClasses;++i)
        {
            // Fuse test is done on average cluster fiber, if it is an active cluster,
            // i.e. at least one of its particle is still moving
            bool activeClass = this->MergeParticleClassFibers(fiberData,mapMergedFibersRef,i);
            if (!activeClass)
                continue;

            for (unsigned int j = i+1;j < numClasses;++j)
            {
                if (classesFusion[j] != j)
                    continue;

                double maxVal = 0;
                bool activeSubClass = this->MergeParticleClassFibers(fiberData,mapMergedFibersFlo,j);
                if (!activeSubClass)
                    continue;

                // Compute a distance between the two clusters, based on user input
                switch (m_ClusterDistance)
                {
                    case 0:
                    {
                        // Former method (quickest)
                        unsigned int minSizeFiber = std::min(mapMergedFibersRef[0].size(),mapMergedFibersFlo[0].size());

                        for (unsigned int l = 0;l < minSizeFiber;++l)
                        {
                            double positionDist = anima::ComputeEuclideanDistance(mapMergedFibersRef[0][l], mapMergedFibersFlo[0][l]);

                            if (positionDist > maxVal)
                                maxVal = positionDist;

                            if (maxVal > m_PositionDistanceFuseThreshold)
                                break;
                        }

                        break;
                    }

                    case 1:
                    {
                        // Hausdorff distance
                        for (unsigned int l = 0;l < mapMergedFibersRef[0].size();++l)
                        {
                            double tmpVal = anima::ComputePointToSetDistance(mapMergedFibersRef[0][l], mapMergedFibersFlo[0]);

                            if (tmpVal > maxVal)
                                maxVal = tmpVal;

                            if (maxVal > m_PositionDistanceFuseThreshold)
                                break;
                        }

                        if (maxVal <= m_PositionDistanceFuseThreshold)
                        {
                            for (unsigned int l = 0;l < mapMergedFibersFlo[0].size();++l)
                            {
                                double tmpVal = anima::ComputePointToSetDistance(mapMergedFibersFlo[0][l], mapMergedFibersRef[0]);

                                if (tmpVal > maxVal)
                                    maxVal = tmpVal;

                                if (maxVal > m_PositionDistanceFuseThreshold)
                                    break;
                            }
                        }

                        break;
                    }

                    case 2:
                    {
                        // Modified Hausdorff distance
                        maxVal = anima::ComputeModifiedDirectedHausdorffDistance(mapMergedFibersRef[0], mapMergedFibersFlo[0]);

                        if (maxVal <= m_PositionDistanceFuseThreshold)
                            maxVal = std::max(maxVal, anima::ComputeModifiedDirectedHausdorffDistance(mapMergedFibersFlo[0], mapMergedFibersRef[0]));

                        break;
                    }

                    default:
                        break;
                }

                // If computed distance is smaller than a threshold, we fuse
                // To do so, an index table (classesFusion) is updated, each of its cells tells
                // to which new class the current class belongs. newNumClasses is the new number of classes
                if (maxVal <= m_PositionDistanceFuseThreshold)
                {
                    newNumClasses--;
                    classesFusion[j] = classesFusion[i];
                }
            }
        }

        // Some post-processing to have contiguous class numbers as an output
        // mapFusion will hold the correspondance between non contiguous and contiguous indexes
        int maxVal = -1;
        unsigned int currentIndex = 0;
        std::map <unsigned int, unsigned int> mapFusion;
        for (unsigned int i = 0;i < numClasses;++i)
        {
            if (maxVal < (int)classesFusion[i])
            {
                mapFusion.insert(std::make_pair(classesFusion[i],currentIndex));
                ++currentIndex;
                maxVal = classesFusion[i];
            }
        }

        for (unsigned int i = 0;i < numClasses;++i)
            classesFusion[i] = mapFusion[classesFusion[i]];
    }

    std::vector <MembershipType> fusedClassesIndexes(newNumClasses);
    for (unsigned int i = 0;i < numClasses;++i)
        fusedClassesIndexes[classesFusion[i]].push_back(i);

    // Now we're done with selecting what to fuse
    // Therefore, deciding on cluster splits. The trick here is that we don't want to actually really perform fuse
    // if it is to split right after (for speed reasons). So we'll play with classesFusion indexes all along
    // to keep track of the original indexes.
    DirectionVectorType afterMergeClassesDirections(newNumClasses);
    MembershipType afterMergeNumPoints(newNumClasses,0);
    std::vector <bool> splitClasses(newNumClasses,false);

    Vector3DType zeroDirection(0.0);
    std::fill(afterMergeClassesDirections.begin(),afterMergeClassesDirections.end(),zeroDirection);

    // Compute average directions after fusion, directions contains the last directions taken by particles
    for (unsigned int i = 0;i < m_NumberOfParticles;++i)
    {
        if (fiberData.particleWeights[i] == 0)
            continue;

        unsigned int classIndex = classesFusion[fiberData.classMemberships[i]];
        for (unsigned int j = 0;j < p;++j)
            afterMergeClassesDirections[classIndex][j] += directions[i][j];
        afterMergeNumPoints[classIndex]++;
    }

    unsigned int numSplits = 0;
    std::vector < std::pair <unsigned int, double> > afterMergeKappaValues;
    // From those averaged directions, we compute a dispersion kappa value, that will be used to decide on split
    for (unsigned int i = 0;i < newNumClasses;++i)
    {
        double norm = 0;
        for (unsigned int j = 0;j < p;++j)
            norm += afterMergeClassesDirections[i][j] * afterMergeClassesDirections[i][j];
        norm = sqrt(norm);

        double R = 1.0;
        double kappa = m_KappaSplitThreshold + 1;
        if (afterMergeNumPoints[i] != 0)
        {
            R = norm / afterMergeNumPoints[i];

            if (R*R > 1.0 - 1.0e-16)
                R = sqrt(1.0 - 1.0e-16);

            kappa = std::exp( anima::safe_log(R) + anima::safe_log(p-R*R) - anima::safe_log(1-R*R) );
        }

        afterMergeKappaValues.push_back(std::make_pair(i,kappa));
        // We do not allow any split resulting in less than m_MinimalNumberOfParticlesPerClass
        // so testing with respect to 2*m_MinimalNumberOfParticlesPerClass
        // If it is ok, and kappa is small enough, there is too much dispersion inside the cluster -> splitting
        if ((kappa <= m_KappaSplitThreshold)&&(afterMergeNumPoints[i] >= 2 * m_MinimalNumberOfParticlesPerClass)&&(afterMergeNumPoints[i] != 0))
            numSplits++;
        else
            afterMergeKappaValues[i].second = m_KappaSplitThreshold + 1;
    }

    std::partial_sort(afterMergeKappaValues.begin(),afterMergeKappaValues.begin() + numSplits,afterMergeKappaValues.end(),pair_comparator());

    for (unsigned int i = 0;i < numSplits;++i)
        splitClasses[afterMergeKappaValues[i].first] = true;

    // Finally apply all this to get our final clusters
    // Each split class will be split into two, so new number of classes
    // after merge and split is newNumClasses + numSplits
    unsigned int finalNumClasses = newNumClasses + numSplits;

    MembershipType newClassesMemberships(m_NumberOfParticles,0);
    std::vector <MembershipType> newReverseClassesMemberships(finalNumClasses);
    MembershipType newClassSizes(finalNumClasses,0);
    ListType newParticleWeights = fiberData.particleWeights;
    ListType newClassWeights(finalNumClasses,0);

    unsigned int currentIndex = 0;

    FiberType vectorToCluster;
    MembershipType clustering;

    // Now, do the real merge/split part
    for (unsigned int i = 0;i < newNumClasses;++i)
    {
        if (!splitClasses[i])
        {
            // ith class is just a potential merge of classes. Easy case: just take all particles
            // from classes marked as new ith class
            for (unsigned int j = 0;j < fusedClassesIndexes[i].size();++j)
            {
                unsigned int classIndex = fusedClassesIndexes[i][j];
                for (unsigned int k = 0;k < fiberData.reverseClassMemberships[classIndex].size();++k)
                {
                    unsigned int particleNumber = fiberData.reverseClassMemberships[classIndex][k];
                    newClassesMemberships[particleNumber] = currentIndex;
                    newReverseClassesMemberships[currentIndex].push_back(particleNumber);
                }
            }

            if (fusedClassesIndexes[i].size() != 1)
            {
                // Recompute class weights after fusion
                newClassWeights[currentIndex] = 0;
                for (unsigned int j = 0;j < fusedClassesIndexes[i].size();++j)
                {
                    unsigned int classIndex = fusedClassesIndexes[i][j];
                    for (unsigned int k = 0;k < fiberData.reverseClassMemberships[classIndex].size();++k)
                        newClassWeights[currentIndex] += fiberData.classWeights[classIndex] * fiberData.particleWeights[fiberData.reverseClassMemberships[classIndex][k]];
                }

                // Recompute particle weights after fusion
                for (unsigned int j = 0;j < fusedClassesIndexes[i].size();++j)
                {
                    unsigned int classIndex = fusedClassesIndexes[i][j];
                    for (unsigned int k = 0;k < fiberData.reverseClassMemberships[classIndex].size();++k)
                    {
                        unsigned int posIndex = fiberData.reverseClassMemberships[classIndex][k];
                        newParticleWeights[posIndex] = std::exp(anima::safe_log(fiberData.classWeights[classIndex]) + anima::safe_log(fiberData.particleWeights[posIndex]) - anima::safe_log(newClassWeights[currentIndex]));
                    }
                }
            }
            else
                newClassWeights[currentIndex] = fiberData.classWeights[fusedClassesIndexes[i][0]];

            ++currentIndex;
        }
        else
        {
            // ith class is a split at the end. In that case, first gather all particles
            // from classes merged before. Then, plug k-means
            vectorToCluster.clear();

            // Gather particles
            for (unsigned int j = 0;j < fusedClassesIndexes[i].size();++j)
            {
                unsigned int classIndex = fusedClassesIndexes[i][j];
                for (unsigned int k = 0;k < fiberData.reverseClassMemberships[classIndex].size();++k)
                    vectorToCluster.push_back(fiberData.fiberParticles[fiberData.reverseClassMemberships[classIndex][k]].back());
            }

            clustering.resize(vectorToCluster.size());

            // Now cluster them, loop until the two classes are not empty
            bool loopOnClustering = true;
            std::uniform_int_distribution <unsigned int> uniInt(0,1);
            while (loopOnClustering)
            {
                for (unsigned int j = 0;j < clustering.size();++j)
                    clustering[j] = uniInt(random_generator) % 2;

                KMeansFilterType kmFilter;
                kmFilter.SetInputData(vectorToCluster);
                kmFilter.SetNumberOfClasses(2);
                kmFilter.InitializeClassesMemberships(clustering);
                kmFilter.SetMaxIterations(100);
                kmFilter.SetVerbose(false);

                kmFilter.Update();

                // Otherwise, go ahead and do the splitting
                clustering = kmFilter.GetClassesMemberships();

                if ((kmFilter.GetNumberPerClass(0) > 0)&&(kmFilter.GetNumberPerClass(1) > 0))
                    loopOnClustering = false;
            }

            // Now assigne new class indexes to particles, plus update class weights
            unsigned int newClassIndex = currentIndex + 1;

            newClassWeights[currentIndex] = 0;
            newClassWeights[newClassIndex] = 0;

            unsigned int pos = 0;
            for (unsigned int j = 0;j < fusedClassesIndexes[i].size();++j)
            {
                unsigned int classIndex = fusedClassesIndexes[i][j];
                for (unsigned int k = 0;k < fiberData.reverseClassMemberships[classIndex].size();++k)
                {
                    unsigned int classPos = currentIndex + clustering[pos];

                    newClassesMemberships[fiberData.reverseClassMemberships[classIndex][k]] = classPos;
                    newReverseClassesMemberships[classPos].push_back(fiberData.reverseClassMemberships[classIndex][k]);

                    newClassWeights[classPos] += fiberData.classWeights[classIndex] * fiberData.particleWeights[fiberData.reverseClassMemberships[classIndex][k]];

                    ++pos;
                }
            }

            // Finally, update particle weights
            pos = 0;
            for (unsigned int j = 0;j < fusedClassesIndexes[i].size();++j)
            {
                unsigned int classIndex = fusedClassesIndexes[i][j];
                for (unsigned int k = 0;k < fiberData.reverseClassMemberships[classIndex].size();++k)
                {
                    unsigned int classPos = currentIndex + clustering[pos];

                    unsigned int posIndex = fiberData.reverseClassMemberships[classIndex][k];
                    newParticleWeights[posIndex] = std::exp(anima::safe_log(fiberData.classWeights[classIndex]) + anima::safe_log(fiberData.particleWeights[posIndex]) - anima::safe_log(newClassWeights[classPos]));

                    ++pos;
                }
            }

            currentIndex += 2;
        }
    }

    double tmpSum = 0;
    for (unsigned int i = 0;i < finalNumClasses;++i)
    {
        if (newReverseClassesMemberships[i].size() > 0)
            newClassWeights[i] = std::max(1.0e-16,newClassWeights[i]);
        tmpSum += newClassWeights[i];
    }

    for (unsigned int i = 0;i < finalNumClasses;++i)
        newClassWeights[i] /= tmpSum;

    for (unsigned int i = 0;i < finalNumClasses;++i)
        newClassSizes[i] = newReverseClassesMemberships[i].size();

    // Replace all fiber data by new ones computed here and we're done
    fiberData.classSizes = newClassSizes;
    fiberData.classWeights = newClassWeights;
    fiberData.classMemberships = newClassesMemberships;
    fiberData.particleWeights = newParticleWeights;
    fiberData.reverseClassMemberships = newReverseClassesMemberships;

    return finalNumClasses;
}

template <class TInputModelImageType>
bool
BaseProbabilisticTractographyImageFilter <TInputModelImageType>
::MergeParticleClassFibers(FiberWorkType &fiberData, FiberProcessVectorType &outputMerged, unsigned int classNumber)
{
    unsigned int numClasses = fiberData.classSizes.size();
    outputMerged.clear();
    if (classNumber >= numClasses)
        return false;

    outputMerged.resize(1);

    std::vector <unsigned int> runningIndexes, stoppedIndexes;
    for (unsigned int j = 0;j < fiberData.classSizes[classNumber];++j)
    {
        if (fiberData.stoppedParticles[fiberData.reverseClassMemberships[classNumber][j]])
            stoppedIndexes.push_back(fiberData.reverseClassMemberships[classNumber][j]);
        else
            runningIndexes.push_back(fiberData.reverseClassMemberships[classNumber][j]);
    }

    FiberType classFiber;
    FiberType tmpFiber;
    unsigned int sizeMerged = 0;
    unsigned int p = PointType::GetPointDimension();

    double sumWeights = 0;
    for (unsigned int j = 0;j < runningIndexes.size();++j)
        sumWeights += fiberData.particleWeights[runningIndexes[j]];

    if (runningIndexes.size() != 0)
    {
        // Use weights provided
        for (unsigned int j = 0;j < runningIndexes.size();++j)
        {
            double tmpWeight = fiberData.particleWeights[runningIndexes[j]];
            if (tmpWeight <= 0)
                continue;

            tmpFiber = fiberData.fiberParticles[runningIndexes[j]];
            for (unsigned int k = 0;k < tmpFiber.size();++k)
            {
                if (k < sizeMerged)
                {
                    for (unsigned int l = 0;l < p;++l)
                        classFiber[k][l] += tmpWeight * tmpFiber[k][l];
                }
                else
                {
                    if (tmpWeight != 0)
                    {
                        sizeMerged++;
                        classFiber.push_back(tmpFiber[k]);
                        for (unsigned int l = 0;l < p;++l)
                            classFiber[k][l] *= tmpWeight;
                    }
                }
            }
        }

        for (unsigned int j = 0;j < sizeMerged;++j)
        {
            for (unsigned int k = 0;k < p;++k)
                classFiber[j][k] /= sumWeights;
        }

        outputMerged[0] = classFiber;
        return true;
    }

    // Treat all fibers equivalently, first construct groups of equal lengths
    std::vector < std::vector <unsigned int> > particleGroups;
    std::vector <unsigned int> particleSizes;
    for (unsigned int i = 0;i < stoppedIndexes.size();++i)
    {
        unsigned int particleSize = fiberData.fiberParticles[stoppedIndexes[i]].size();
        bool sizeFound = false;
        for (unsigned int j = 0;j < particleSizes.size();++j)
        {
            if (particleSize == particleSizes[j])
            {
                particleGroups[j].push_back(stoppedIndexes[i]);
                sizeFound = true;
                break;
            }
        }

        if (!sizeFound)
        {
            particleSizes.push_back(particleSize);
            std::vector <unsigned int> tmpVec(1,stoppedIndexes[i]);
            particleGroups.push_back(tmpVec);
        }
    }

    // For each group of equal length, build fiber
    outputMerged.resize(particleGroups.size());
    for (unsigned int i = 0;i < particleGroups.size();++i)
    {
        classFiber.clear();
        sizeMerged = 0;

        for (unsigned int j = 0;j < particleGroups[i].size();++j)
        {
            tmpFiber = fiberData.fiberParticles[particleGroups[i][j]];
            for (unsigned int k = 0;k < tmpFiber.size();++k)
            {
                if (k < sizeMerged)
                {
                    for (unsigned int l = 0;l < p;++l)
                        classFiber[k][l] += tmpFiber[k][l];
                }
                else
                {
                    sizeMerged++;
                    classFiber.push_back(tmpFiber[k]);
                }
            }
        }

        for (unsigned int j = 0;j < sizeMerged;++j)
        {
            for (unsigned int k = 0;k < p;++k)
                classFiber[j][k] /= particleGroups[i].size();
        }

        outputMerged[i] = classFiber;
    }

    return false;
}

} // end of namespace anima
