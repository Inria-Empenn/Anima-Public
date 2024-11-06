#pragma once
#include "animaMultiThreadedSTAPLEImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkPoolMultiThreader.h>
#include <itkTimeProbe.h>

namespace anima
{

template <typename TInputImage>
void
MultiThreadedSTAPLEImageFilter <TInputImage>
::CheckComputationMask()
{
    if (this->GetComputationMask())
        return;

    MaskImagePointer tmpMask = MaskImageType::New();
    tmpMask->Initialize();

    tmpMask->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
    tmpMask->SetSpacing(this->GetInput(0)->GetSpacing());
    tmpMask->SetOrigin(this->GetInput(0)->GetOrigin());
    tmpMask->SetDirection(this->GetInput(0)->GetDirection());

    tmpMask->Allocate();
    tmpMask->FillBuffer(1);
    this->SetComputationMask(tmpMask);

    if (m_MaskDilationRadius == 0)
        return;

    this->GetComputationMask()->FillBuffer(0);

    typedef itk::ImageRegionConstIterator <TInputImage> InIteratorType;
    typedef itk::ImageRegionIterator <MaskImageType> MaskRegionIteratorType;

    MaskRegionIteratorType maskItr(this->GetComputationMask(),this->GetInput(0)->GetLargestPossibleRegion());
    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
    {
        InIteratorType itr(this->GetInput(i),this->GetInput(i)->GetLargestPossibleRegion());

        maskItr.GoToBegin();
        while (!maskItr.IsAtEnd())
        {
            if (maskItr.Get() != 0)
            {
                ++itr;
                ++maskItr;
                continue;
            }

            if (itr.Get() != 0)
                maskItr.Set(1);

            ++itr;
            ++maskItr;
        }
    }

    typedef itk::BinaryBallStructuringElement <unsigned short, 3> BallElementType;
    typedef itk::GrayscaleDilateImageFilter <MaskImageType,MaskImageType,BallElementType> DilateFilterType;

    typename DilateFilterType::Pointer dilateFilter = DilateFilterType::New();
    dilateFilter->SetInput(this->GetComputationMask());
    dilateFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

    BallElementType tmpBall;
    BallElementType::SizeType ballSize;
    ballSize[0] = m_MaskDilationRadius;
    ballSize[1] = m_MaskDilationRadius;
    ballSize[2] = m_MaskDilationRadius;
    tmpBall.SetRadius(ballSize);
    tmpBall.CreateStructuringElement();

    dilateFilter->SetKernel(tmpBall);
    dilateFilter->Update();

    typename MaskImageType::Pointer dilMask = dilateFilter->GetOutput();
    dilMask->DisconnectPipeline();

    this->SetComputationMask(dilMask);
}

template <typename TInputImage>
void
MultiThreadedSTAPLEImageFilter <TInputImage>
::InitializeMissingStructures()
{
    if (m_nbClasses <= 0)
        this->InitializeNbClassesFromData();

    m_MAPStructures.clear();
    std::vector <unsigned int> tmpVec(m_nbClasses,1);
    if (!m_AccountForMissingStructures)
    {
        for (unsigned int i = 0;i < this->GetNumberOfInputs();++i)
            m_MAPStructures.push_back(tmpVec);

        return;
    }

    typedef itk::ImageRegionConstIterator <TInputImage> InIteratorType;
    typedef itk::ImageRegionIterator <MaskImageType> MaskRegionIteratorType;

    InputImageRegionType largestRegion = this->GetInput(0)->GetLargestPossibleRegion();
    MaskRegionIteratorType maskItr(this->GetComputationMask(),largestRegion);

    for (unsigned int i = 0;i < this->GetNumberOfInputs();++i)
    {
        for (unsigned int j = 0;j < m_nbClasses;++j)
            tmpVec[j] = 0;

        InIteratorType tmpIt(this->GetInput(i), largestRegion);
        maskItr.GoToBegin();

        unsigned int sumStrs = 0;
        while (!maskItr.IsAtEnd())
        {
            if (maskItr.Get() != 0)
            {
                unsigned int tmpVal = tmpIt.Get();
                if (tmpVec[tmpVal] == 0)
                {
                    sumStrs++;
                    tmpVec[tmpVal] = 1;
                }
            }

            if (sumStrs == m_nbClasses)
                break;

            ++maskItr;
            ++tmpIt;
        }

        m_MAPStructures.push_back(tmpVec);
    }
}

template <typename TInputImage>
void
MultiThreadedSTAPLEImageFilter <TInputImage>
::InitializeExpertParameters(double diagValue)
{
    this->CheckComputationMask();

    if (m_nbClasses <= 0)
        this->InitializeNbClassesFromData();

    if (m_MAPStructures.size() == 0)
        this->InitializeMissingStructures();

    // m_nbClasses has to have been computed from the images
    unsigned int nbExperts = this->GetNumberOfIndexedInputs();

    double nonDiagValue = (1.0 - diagValue)/(m_nbClasses - 1);

    m_OldExpParams.clear();
    m_ExpParams.clear();

    ParametersMatrixType tmpMat(m_nbClasses,m_nbClasses);
    for (unsigned int i = 0;i < nbExperts;++i)
    {
        for (unsigned int k = 0;k < m_nbClasses;++k)
        {
            if (m_MAPStructures[i][k] != 0)
            {
                for (unsigned int j = 0;j < m_nbClasses;++j)
                {
                    if (j != k)
                        tmpMat(j,k) = nonDiagValue;
                    else
                        tmpMat(j,j) = diagValue;
                }
            }
            else
            {
                // Missing structure
                tmpMat(0,k) = diagValue;
                for (unsigned int j = 1;j < m_nbClasses;++j)
                    tmpMat(j,k) = nonDiagValue;
            }
        }

        m_OldExpParams.push_back(tmpMat);
        m_ExpParams.push_back(tmpMat);
    }
}

template <typename TInputImage>
bool
MultiThreadedSTAPLEImageFilter <TInputImage>
::endConditionReached()
{
    double absDiff = 0;

    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
    {
        for (unsigned int j = 0;j < m_nbClasses;++j)
            for (unsigned int k = 0;k < m_nbClasses;++k)
            {
                if (fabs(m_OldExpParams[i](j,k) - m_ExpParams[i](j,k)) > absDiff)
                    absDiff = fabs(m_OldExpParams[i](j,k) - m_ExpParams[i](j,k));
            }
    }

    if (absDiff > m_RelativeConvergenceThreshold)
        return false;
    else
        return true;
}

template <typename TInputImage>
void
MultiThreadedSTAPLEImageFilter <TInputImage>
::InitializeNbClassesFromData()
{
    // In the future implement here a way to rescale non contiguous labels so that they are after
    // Right now just get the max value over all images and set it as nbClasses
    this->CheckComputationMask();

    m_nbClasses = 0;

    typedef itk::ImageRegionConstIteratorWithIndex <TInputImage> InIteratorType;
    typedef itk::ImageRegionIteratorWithIndex <MaskImageType> MaskRegionIteratorType;

    InputImageRegionType largestRegion = this->GetInput(0)->GetLargestPossibleRegion();
    MaskRegionIteratorType maskItr(this->GetComputationMask(),largestRegion);

    for (unsigned int i = 0;i < this->GetNumberOfInputs();++i)
    {
        maskItr.GoToBegin();
        InIteratorType inItr(this->GetInput(i),largestRegion);
        inItr.GoToBegin();

        while(!maskItr.IsAtEnd())
        {
            if (maskItr.Get() != 0)
            {
                unsigned int tmpVal = inItr.Get();
                m_nbClasses = std::max(m_nbClasses,tmpVal);
            }

            ++maskItr;
            ++inItr;
        }
    }

    if (m_nbClasses != 0)
        m_nbClasses++;
}

template <typename TInputImage>
void
MultiThreadedSTAPLEImageFilter <TInputImage>
::InitializePriorFromData()
{
    this->CheckComputationMask();

    if (m_nbClasses <= 0)
        this->InitializeNbClassesFromData();

    if (m_MAPStructures.size() != this->GetNumberOfInputs())
        this->InitializeMissingStructures();

    if (!m_GTPriorImage.IsNull())
    {
        if (m_GTPriorImage->GetNumberOfComponentsPerPixel() == m_nbClasses)
            return;
        else
            itkExceptionMacro("Number of classes in ground truth prior image is not the same as in the images...");
    }

    m_Prior.resize(m_nbClasses);

    for (unsigned int i = 0;i < m_nbClasses;++i)
        m_Prior[i] = 0;

    typedef itk::ImageRegionConstIteratorWithIndex <TInputImage> InIteratorType;
    typedef itk::ImageRegionIteratorWithIndex <MaskImageType> MaskRegionIteratorType;

    InputImageRegionType largestRegion = this->GetInput(0)->GetLargestPossibleRegion();

    MaskRegionIteratorType maskItr(this->GetComputationMask(),largestRegion);
    double nbRealPts = 0;

    for (unsigned int i = 0;i < this->GetNumberOfInputs();++i)
    {
        maskItr.GoToBegin();
        InIteratorType inputIterator(this->GetInput(i), largestRegion);
        inputIterator.GoToBegin();

        while (!maskItr.IsAtEnd())
        {
            if (maskItr.Get() != 0)
            {
                m_Prior[inputIterator.Get()]++;
                if (i == 0)
                    nbRealPts++;
            }

            ++maskItr;
            ++inputIterator;
        }
    }

    double sumPriors = 0;

    for (unsigned int i = 1;i < m_nbClasses;++i)
    {
        double nbExpClass = 0;
        for (unsigned int j = 0;j < this->GetNumberOfInputs();++j)
            nbExpClass += (m_MAPStructures[j][i] != 0);
        m_Prior[i] /= (nbRealPts*nbExpClass);
        sumPriors += m_Prior[i];
    }

    m_Prior[0] = 1 - sumPriors;
}

template <typename TInputImage>
void
MultiThreadedSTAPLEImageFilter <TInputImage>
::GenerateOutputInformation()
{
    // Override the method in itkImageSource, so we can set the vector length of
    // the output itk::VectorImage

    this->Superclass::GenerateOutputInformation();

    if (m_nbClasses <= 0)
        this->InitializeNbClassesFromData();

    OutputImageType *output = this->GetOutput();
    output->SetVectorLength(m_nbClasses);
}

template <typename TInputImage>
void
MultiThreadedSTAPLEImageFilter <TInputImage>
::GenerateData()
{
    this->AllocateOutputs();
    this->BeforeThreadedGenerateData();

    if (m_MAPUpdate == DIAGONAL_MAP)
        this->SetAccountForMissingStructures(false);

    this->InitializeMissingStructures();

    if (m_nbClasses != m_Prior.size())
        this->InitializePriorFromData();

    if ((m_nbClasses == 2)&&(m_MAPUpdate == FULL_MAP))
    {
        m_AlphaMAP += m_BetaMAPNonDiag - 1;
        m_BetaMAP += m_AlphaMAPNonDiag - 1;
        m_MAPUpdate = DIAGONAL_MAP;
    }

    if ((m_MAPUpdate == FULL_MAP)&&(m_AlphaMAPNonDiag == 1)&&(m_BetaMAPNonDiag == 1))
        m_MAPUpdate = DIAGONAL_MAP;

    if (this->GetNumberOfInputs() != m_ExpParams.size())
        this->InitializeExpertParameters(0.99);

    unsigned int itncount = 0;
    bool continueLoop = true;
    while ((itncount < m_MaximumIterations)&&(continueLoop))
    {
        if (m_Verbose)
            std::cout << "Iteration " << itncount + 1 << "..." << std::endl;

        itk::TimeProbe tmpTime;
        tmpTime.Start();

        // Parallelize by calling DynamicThreadedGenerateData, estimation of reference standard
        this->GetMultiThreader()->template ParallelizeImageRegion<TOutputImage::ImageDimension> (
            this->GetOutput()->GetRequestedRegion(),
            [this](const OutputImageRegionType & outputRegionForThread)
              { this->DynamicThreadedGenerateData(outputRegionForThread); }, this);

        tmpTime.Stop();

        if (m_Verbose)
            std::cout << "Reference standard estimated in " << tmpTime.GetTotal() << "..." << std::endl;

        itk::TimeProbe tmpTimePerf;
        tmpTimePerf.Start();

        EstimatePerformanceParameters();

        tmpTimePerf.Stop();

        if (m_Verbose)
            std::cout << "Performance parameters estimated in " << tmpTimePerf.GetTotal() << "..." << std::endl;

        ++itncount;

        if (itncount != 1)
            continueLoop = !endConditionReached();

        if (continueLoop)
        {
            for (unsigned int i = 0;i < this->GetNumberOfInputs();++i)
            {
                for (unsigned int j = 0;j < m_nbClasses;++j)
                    for (unsigned int k = 0;k < m_nbClasses;++k)
                        m_OldExpParams[i](j,k) = m_ExpParams[i](j,k);
            }
        }
    }

    m_ElapsedIterations = itncount;
}

template <typename TInputImage>
void
MultiThreadedSTAPLEImageFilter <TInputImage>
::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    unsigned int nbExperts = this->GetNumberOfInputs();

    typedef itk::ImageRegionConstIteratorWithIndex <TInputImage> InIteratorType;
    typedef itk::ImageRegionIteratorWithIndex <TOutputImage> OutRegionIteratorType;
    typedef itk::ImageRegionConstIteratorWithIndex <GTPriorImageType> GTPriorIteratorType;
    typedef itk::ImageRegionIteratorWithIndex <MaskImageType> MaskRegionIteratorType;

    OutRegionIteratorType outItr(this->GetOutput(),outputRegionForThread);

    GTPriorIteratorType gtPriorIt;
    GTPriorImageType::PixelType lPrior(m_nbClasses);

    if (m_GTPriorImage.IsNull())
    {
        for (unsigned int i = 0;i < m_nbClasses;++i)
            lPrior[i] = m_Prior[i];
    }
    else
    {
        gtPriorIt = GTPriorIteratorType(m_GTPriorImage,outputRegionForThread);
    }

    std::vector <InIteratorType> inputIterators;
    for (unsigned int i = 0; i < nbExperts;++i)
    {
        inputIterators.push_back(InIteratorType(this->GetInput(i), outputRegionForThread));
        inputIterators[i].GoToBegin();
    }

    MaskRegionIteratorType maskItr(this->GetComputationMask(),outputRegionForThread);

    outItr.GoToBegin();
    maskItr.GoToBegin();

    OutputPixelType tmpClassif(m_nbClasses);

    while (!outItr.IsAtEnd())
    {
        if (maskItr.Get() == 0)
        {
            ++outItr;
            for (unsigned int i = 0; i < nbExperts;++i)
                ++inputIterators[i];
            ++maskItr;

            if (!m_GTPriorImage.IsNull())
                ++gtPriorIt;

            continue;
        }

        if (!m_GTPriorImage.IsNull())
            lPrior = gtPriorIt.Get();

        long double denom = 0;
        for (unsigned int m = 0;m < m_nbClasses;++m)
        {
            long double temp = 1;
            for (unsigned int k = 0;k < nbExperts;++k)
            {
                unsigned int l = inputIterators[k].Get();
                temp *= m_ExpParams[k](l,m);
            }
            denom += lPrior[m]*temp;
        }

        for (unsigned int m = 0;m < m_nbClasses;++m)
        {
            long double pDijTi = 1;
            for (unsigned int k = 0;k < nbExperts;++k)
            {
                unsigned int l = inputIterators[k].Get();
                pDijTi *= m_ExpParams[k](l,m);
            }
            pDijTi *= lPrior[m];

            tmpClassif[m] = pDijTi/denom;
        }

        outItr.Set(tmpClassif);

        ++outItr;
        ++maskItr;
        for (unsigned int i = 0; i < nbExperts;++i)
            ++inputIterators[i];

        if (!m_GTPriorImage.IsNull())
            ++gtPriorIt;
    }
}

template <typename TInputImage>
void
MultiThreadedSTAPLEImageFilter <TInputImage>
::EstimatePerformanceParameters()
{
    itk::PoolMultiThreader::Pointer threaderMstep = itk::PoolMultiThreader::New();

    EMStepThreadStruct *tmpStr = new EMStepThreadStruct;
    tmpStr->Filter = this;

    unsigned int actualNumberOfThreads = std::min(this->GetNumberOfWorkUnits(),(unsigned int)this->GetNumberOfInputs());

    threaderMstep->SetNumberOfWorkUnits(actualNumberOfThreads);
    threaderMstep->SetSingleMethod(this->ThreadEstimatePerfParams,tmpStr);
    threaderMstep->SingleMethodExecute();

    delete tmpStr;
}

template <typename TInputImage>
ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION
MultiThreadedSTAPLEImageFilter <TInputImage>
::ThreadEstimatePerfParams(void *arg)
{
    itk::MultiThreaderBase::WorkUnitInfo *threadArgs = (itk::MultiThreaderBase::WorkUnitInfo *)arg;

    unsigned int nbThread = threadArgs->WorkUnitID;
    unsigned int nbProcs = threadArgs->NumberOfWorkUnits;

    EMStepThreadStruct *tmpStr = (EMStepThreadStruct *)threadArgs->UserData;
    unsigned int nbExperts = tmpStr->Filter->GetNumberOfInputs();

    unsigned int minExp = (unsigned int)floor((double)nbThread*nbExperts/nbProcs);
    unsigned int maxExp = (unsigned int)floor((double)(nbThread + 1.0)*nbExperts/nbProcs);

    maxExp = std::min(nbExperts,maxExp);

    switch (tmpStr->Filter->GetMAPUpdate())
    {
        case FULL_MAP:
            tmpStr->Filter->EstimateFullMAPPerformanceParameters(minExp,maxExp);
            break;

        case DIAGONAL_MAP:
            tmpStr->Filter->EstimateMAPPerformanceParameters(minExp,maxExp);
            break;

        case STANDARD:
        default:
            tmpStr->Filter->EstimatePerformanceParameters(minExp,maxExp);
            break;
    }

    return ITK_THREAD_RETURN_DEFAULT_VALUE;
}

template <typename TInputImage>
void
MultiThreadedSTAPLEImageFilter <TInputImage>
::EstimatePerformanceParameters(unsigned int minExp, unsigned int maxExp)
{
    typedef itk::ImageRegionConstIteratorWithIndex <TInputImage> InIteratorType;
    typedef itk::ImageRegionIteratorWithIndex <TOutputImage> OutRegionIteratorType;
    typedef itk::ImageRegionIteratorWithIndex <MaskImageType> MaskRegionIteratorType;

    OutRegionIteratorType outItr(this->GetOutput(),this->GetOutput()->GetRequestedRegion());
    MaskRegionIteratorType maskItr(this->GetComputationMask(),this->GetOutput()->GetRequestedRegion());

    for (unsigned int i = minExp; i < maxExp;++i)
    {
        InIteratorType inputIterator(this->GetInput(i),this->GetOutput()->GetRequestedRegion());
        std::vector <std::vector <long double> > nums;
        std::vector <long double> denom(m_nbClasses,0);

        inputIterator.GoToBegin();
        outItr.GoToBegin();
        maskItr.GoToBegin();

        for (unsigned int j = 0;j < m_nbClasses;++j)
            nums.push_back(std::vector<long double> (m_nbClasses,0));

        while (!inputIterator.IsAtEnd())
        {
            if (maskItr.Get() != 0)
            {
                for (unsigned int j = 0;j < m_nbClasses;++j)
                    nums[j][inputIterator.Get()] += outItr.Get()[j];
            }

            ++maskItr;
            ++outItr;
            ++inputIterator;
        }

        for (unsigned int j = 0;j < m_nbClasses;++j)
            for (unsigned int k = 0;k < m_nbClasses;++k)
                denom[j] += nums[j][k];

        for (unsigned int j = 0;j < m_nbClasses;++j)
            for (unsigned int k = 0;k < m_nbClasses;++k)
                m_ExpParams[i](k,j) = nums[j][k]/denom[j];
    }
}

template <typename TInputImage>
void
MultiThreadedSTAPLEImageFilter <TInputImage>
::EstimateMAPPerformanceParameters(unsigned int minExp, unsigned int maxExp)
{
    typedef itk::ImageRegionConstIteratorWithIndex <TInputImage> InIteratorType;
    typedef itk::ImageRegionIteratorWithIndex <TOutputImage> OutRegionIteratorType;
    typedef itk::ImageRegionIteratorWithIndex <MaskImageType> MaskRegionIteratorType;

    OutRegionIteratorType outItr(this->GetOutput(),this->GetOutput()->GetRequestedRegion());
    MaskRegionIteratorType maskItr(this->GetComputationMask(),this->GetOutput()->GetRequestedRegion());

    for (unsigned int i = minExp; i < maxExp;++i)
    {
        InIteratorType inputIterator(this->GetInput(i),this->GetOutput()->GetRequestedRegion());
        std::vector <std::vector <long double> > nums;
        std::vector <long double> denom(m_nbClasses,0);

        inputIterator.GoToBegin();
        outItr.GoToBegin();
        maskItr.GoToBegin();

        for (unsigned int j = 0;j < m_nbClasses;++j)
            nums.push_back(std::vector<long double> (m_nbClasses,0));

        while (!inputIterator.IsAtEnd())
        {
            if (maskItr.Get() != 0)
            {
                for (unsigned int j = 0;j < m_nbClasses;++j)
                    nums[j][inputIterator.Get()] += outItr.Get()[j];
            }

            ++maskItr;
            ++outItr;
            ++inputIterator;
        }

        for (unsigned int j = 0;j < m_nbClasses;++j)
            for (unsigned int k = 0;k < m_nbClasses;++k)
                denom[j] += nums[j][k];

        for (unsigned int j = 0;j < m_nbClasses;++j)
            for (unsigned int k = 0;k < m_nbClasses;++k)
            {
                if (k != j)
                {
                    if (denom[j] != nums[j][j])
                    {
                        m_ExpParams[i](k,j) = nums[j][k]*(m_MAPWeighting*(m_BetaMAP - 1) + denom[j] - nums[j][j])/
                                ((denom[j] + m_MAPWeighting*(m_AlphaMAP + m_BetaMAP - 2))*(denom[j] - nums[j][j]));
                    }
                    else
                    {
                        if (m_nbClasses == 2)
                            m_ExpParams[i](k,j) = (m_MAPWeighting*(m_BetaMAP - 1) + denom[j] - nums[j][j])/
                                    (denom[j] + m_MAPWeighting*(m_AlphaMAP + m_BetaMAP - 2));
                        else
                            m_ExpParams[i](k,j) = 0;
                    }
                }
                else
                    m_ExpParams[i](k,j) = (nums[j][k] + m_MAPWeighting*(m_AlphaMAP - 1))/
                            (denom[j] + m_MAPWeighting*(m_AlphaMAP + m_BetaMAP - 2));
            }
    }
}

template <typename TInputImage>
void
MultiThreadedSTAPLEImageFilter <TInputImage>
::EstimateFullMAPPerformanceParameters(unsigned int minExp, unsigned int maxExp)
{
    typedef itk::ImageRegionConstIteratorWithIndex <TInputImage> InIteratorType;
    typedef itk::ImageRegionIteratorWithIndex <TOutputImage> OutRegionIteratorType;
    typedef itk::ImageRegionIteratorWithIndex <MaskImageType> MaskRegionIteratorType;

    OutRegionIteratorType outItr(this->GetOutput(),this->GetOutput()->GetRequestedRegion());
    MaskRegionIteratorType maskItr(this->GetComputationMask(),this->GetOutput()->GetRequestedRegion());

    long double mapDenom = m_MAPWeighting*((m_nbClasses - 1)*(m_AlphaMAPNonDiag + m_BetaMAPNonDiag - 2) + m_AlphaMAP + m_BetaMAP - 2);

    for (unsigned int i = minExp; i < maxExp;++i)
    {
        InIteratorType inputIterator(this->GetInput(i),this->GetOutput()->GetRequestedRegion());
        std::vector <std::vector <long double> > nums;
        std::vector <long double> denom(m_nbClasses,mapDenom);

        inputIterator.GoToBegin();
        outItr.GoToBegin();
        maskItr.GoToBegin();

        for (unsigned int j = 0;j < m_nbClasses;++j)
            nums.push_back(std::vector<long double> (m_nbClasses,0));

        while (!inputIterator.IsAtEnd())
        {
            if (maskItr.Get() != 0)
            {
                for (unsigned int j = 0;j < m_nbClasses;++j)
                    nums[j][inputIterator.Get()] += outItr.Get()[j];
            }

            ++maskItr;
            ++outItr;
            ++inputIterator;
        }

        for (unsigned int j = 0;j < m_nbClasses;++j)
            for (unsigned int k = 0;k < m_nbClasses;++k)
                denom[j] += nums[j][k];

        for (unsigned int j = 0;j < m_nbClasses;++j)
            this->ComputeFixedPointMAPEstimates(i,j,nums[j],denom[j]);
    }
}

template <typename TInputImage>
void
MultiThreadedSTAPLEImageFilter <TInputImage>
::ComputeFixedPointMAPEstimates(unsigned int numExp, unsigned int numClass, std::vector <long double> &constantNums, long double &constantDenom)
{
    std::vector <long double> newEstimates(m_nbClasses,0), oldEstimates(m_nbClasses,0);

    for (unsigned int i = 0;i < m_nbClasses;++i)
        oldEstimates[i] = m_ExpParams[numExp](i,numClass);

    unsigned int curIt = 0;
    long double quickDiff = 10.0;

    while (((curIt < m_MaxIterFixedPoint)&&(quickDiff > m_EpsilonFixedPoint))||(curIt == 0))
    {
        curIt++;

        for (unsigned int lp = 0;lp < m_nbClasses;++lp)
        {
            long double variablePart = 0;
            for (unsigned int np = 0;np < m_nbClasses;++np)
            {
                if (np != lp)
                {
                    if (m_MAPStructures[numExp][numClass] == 1)
                    {
                        if (np != numClass)
                            variablePart += (m_BetaMAPNonDiag - 1) / (oldEstimates[np] - 1);
                        else
                            variablePart += (m_BetaMAP - 1) / (oldEstimates[np] - 1);
                    }
                    else
                    {
                        // Missing structure case, the diagonal term is transferred on theta_{j0l}
                        if (np != 0)
                            variablePart += (m_BetaMAPNonDiag - 1) / (oldEstimates[np] - 1);
                        else
                            variablePart += (m_BetaMAP - 1) / (oldEstimates[np] - 1);
                    }
                }
            }

            if (m_MAPStructures[numExp][numClass] == 1)
            {
                if (lp != numClass)
                    newEstimates[lp] = constantNums[lp] + m_MAPWeighting*(m_AlphaMAPNonDiag - 1);
                else
                    newEstimates[lp] = constantNums[lp] + m_MAPWeighting*(m_AlphaMAP - 1);
            }
            else
            {
                // Missing structure case, the diagonal term is transferred on theta_{j0l}
                if (lp != 0)
                    newEstimates[lp] = constantNums[lp] + m_MAPWeighting*(m_AlphaMAPNonDiag - 1);
                else
                    newEstimates[lp] = constantNums[lp] + m_MAPWeighting*(m_AlphaMAP - 1);
            }

            newEstimates[lp] /= (constantDenom + m_MAPWeighting*variablePart);
        }

        quickDiff = 0;

        for (unsigned int lp = 0;lp < m_nbClasses;++lp)
            quickDiff = std::max(quickDiff,std::abs(newEstimates[lp] - oldEstimates[lp]));

        if (quickDiff > m_EpsilonFixedPoint)
        {
            for (unsigned int i = 0;i < m_nbClasses;++i)
                oldEstimates[i] = newEstimates[i];
        }
    }

    for (unsigned int i = 0;i < m_nbClasses;++i)
        m_ExpParams[numExp](i,numClass) = newEstimates[i];
}

template <typename TInputImage>
void
MultiThreadedSTAPLEImageFilter <TInputImage>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
    os << indent << "Expert parameters:" << std::endl;
    for (unsigned int i = 0;i < this->GetNumberOfInputs();++i)
    {
        os << indent << "Expert " << i << " : " << std::endl;
        for (unsigned int j = 0;j < m_nbClasses;++j)
        {
            for (unsigned int k = 0;k < m_nbClasses;++k)
                os << indent << m_ExpParams[i](j,k) << " ";

            os << indent << std::endl;
        }

        os << indent << std::endl;
    }

    os << indent << std::endl;

    if (m_GTPriorImage.IsNull())
    {
        os << indent << "Prior values:" << std::endl;

        for (unsigned int i = 0;i < m_Prior.size();++i)
            os << indent << m_Prior[i] << " ";

        os << indent << std::endl;
    }
}

template <typename TInputImage>
TInputImage *
MultiThreadedSTAPLEImageFilter <TInputImage>
::GetClassificationAsLabelMap()
{
    if (m_LabelMap)
        m_LabelMap->Delete();

    m_LabelMap = InputImageType::New();

    m_LabelMap->Initialize();
    m_LabelMap->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
    m_LabelMap->SetSpacing(this->GetInput(0)->GetSpacing());
    m_LabelMap->SetOrigin(this->GetInput(0)->GetOrigin());
    m_LabelMap->SetDirection(this->GetInput(0)->GetDirection());

    m_LabelMap->Allocate();

    typedef itk::ImageRegionIteratorWithIndex <TInputImage> InIteratorType;
    typedef itk::ImageRegionIteratorWithIndex <TOutputImage> OutRegionIteratorType;
    typedef itk::ImageRegionIteratorWithIndex <MaskImageType> MaskRegionIteratorType;

    InputImageRegionType largestRegion = this->GetOutput()->GetLargestPossibleRegion();
    InIteratorType labelMapItr(m_LabelMap,largestRegion);
    OutRegionIteratorType outItr(this->GetOutput(),largestRegion);

    MaskRegionIteratorType maskItr(this->GetComputationMask(),largestRegion);

    while(!maskItr.IsAtEnd())
    {
        if (maskItr.Get() != 0)
        {
            OutputPixelType tmpVal = outItr.Get();
            unsigned int maxClass = 0;
            double maxVal = tmpVal[0];

            for (unsigned int i = 0;i < m_nbClasses;++i)
            {
                if (tmpVal[i] > maxVal)
                {
                    maxVal = tmpVal[i];
                    maxClass = i;
                }
            }

            labelMapItr.Set(maxClass);
        }
        else
            labelMapItr.Set(0);

        ++labelMapItr;
        ++outItr;
        ++maskItr;
    }

    return m_LabelMap;
}

template <typename TInputImage>
itk::Image <double, 3> *
MultiThreadedSTAPLEImageFilter <TInputImage>
::GetExpertParametersAsImage()
{
    ParamsImageType::Pointer resVal = ParamsImageType::New();
    InputImageRegionType region;
    region.SetIndex(0,0);
    region.SetIndex(1,0);
    region.SetIndex(2,0);

    region.SetSize(0,this->GetNumberOfInputs());
    region.SetSize(1,m_nbClasses);
    region.SetSize(2,m_nbClasses);

    resVal->Initialize();

    resVal->SetRegions (region);
    resVal->Allocate();

    InputImageIndexType tmpInd;
    for (unsigned int i = 0;i < this->GetNumberOfInputs();++i)
    {
        tmpInd[0] = i;
        for (unsigned j = 0;j < m_nbClasses;++j)
        {
            tmpInd[1] = j;
            for (unsigned k = 0;k < m_nbClasses;++k)
            {
                tmpInd[2] = k;
                resVal->SetPixel(tmpInd,m_ExpParams[i](j,k));
            }
        }
    }

    resVal->Register();

    return resVal;
}

} // end namespace anima
