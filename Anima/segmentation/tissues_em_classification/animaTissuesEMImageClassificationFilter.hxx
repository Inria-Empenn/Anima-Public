#pragma once
#include "animaTissuesEMClassificationImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkPoolMultiThreader.h>
#include <itkTimeProbe.h>

#include <animaBaseTensorTools.h>

namespace anima
{

template <typename TInputImage>
bool
TissuesEMClassificationImageFilter <TInputImage>
::endConditionReached()
{
    for (unsigned int i = 0;i < m_NumberOfClasses;++i)
    {
        for (unsigned int j = 0;j < m_NumberOfInputs;++j)
        {
            double denom = (std::abs(m_OldClassesMeans[i][j]) + std::abs(m_ClassesMeans[i][j])) / 2.0;
            if (denom > 0)
            {
                double absDiff = std::abs(m_OldClassesMeans[i][j] - m_ClassesMeans[i][j]);
                if (absDiff / denom > m_RelativeConvergenceThreshold)
                    return false;
            }
        }

        for (unsigned int j = 0;j < m_NumberOfInputs;++j)
        {
            for (unsigned int k = j;k < m_NumberOfInputs;++k)
            {
                double denom = (std::abs(m_OldClassesVariances[i](j,k)) + std::abs(m_ClassesVariances[i](j,k))) / 2.0;
                if (denom > 0)
                {
                    double absDiff = std::abs(m_OldClassesVariances[i](j,k) - m_ClassesVariances[i](j,k));
                    if (absDiff / denom > m_RelativeConvergenceThreshold)
                        return false;
                }
            }
        }
    }

    return true;
}

template <typename TInputImage>
void
TissuesEMClassificationImageFilter <TInputImage>
::GenerateOutputInformation()
{
    // Override the method in itkImageSource, so we can set the vector length of
    // the output itk::VectorImage

    this->Superclass::GenerateOutputInformation();

    m_NumberOfClasses = m_LocalPriorImage->GetNumberOfComponentsPerPixel();
    m_NumberOfInputs = this->GetNumberOfIndexedInputs();

    OutputImageType *output = this->GetOutput();
    output->SetVectorLength(m_NumberOfClasses);
}

template <typename TInputImage>
void
TissuesEMClassificationImageFilter <TInputImage>
::BeforeThreadedGenerateData()
{
    this->Superclass::BeforeThreadedGenerateData();

    if (m_LocalPriorImage.IsNull())
        itkExceptionMacro("A local prior image is necessary for this implementation");

    typedef itk::ImageRegionIterator <LocalPriorImageType> LocalPriorIteratorType;
    typedef itk::ImageRegionIterator <OutputImageType> OutIteratorType;
    typedef itk::ImageRegionIterator <MaskImageType> MaskRegionIteratorType;

    LocalPriorIteratorType localPriorItr(m_LocalPriorImage,this->GetComputationRegion());
    OutIteratorType outItr(this->GetOutput(),this->GetComputationRegion());
    MaskRegionIteratorType maskItr(this->GetComputationMask(),this->GetComputationRegion());

    OutputPixelType outPixel(m_NumberOfClasses);

    while(!maskItr.IsAtEnd())
    {
        if (maskItr.Get() != 0)
        {
            outPixel = localPriorItr.Get();
            double sumProbas = 0.0;
            for (unsigned int i = 0;i < m_NumberOfClasses;++i)
                sumProbas += outPixel[i];

            if (sumProbas == 0)
            {
                outPixel.Fill(1.0 / m_NumberOfClasses);
                localPriorItr.Set(outPixel);
            }

            outItr.Set(outPixel);
        }

        ++maskItr;
        ++localPriorItr;
        ++outItr;
    }

    m_ClassesMeans.resize(m_NumberOfClasses);
    m_OldClassesMeans.resize(m_NumberOfClasses);
    m_ClassesVariances.resize(m_NumberOfClasses);
    m_InverseClassesVariances.resize(m_NumberOfClasses);
    m_ClassesVariancesSqrtDeterminants.resize(m_NumberOfClasses);
    m_OldClassesVariances.resize(m_NumberOfClasses);

    for (unsigned int i = 0;i < m_NumberOfClasses;++i)
    {
        m_ClassesVariances[i].set_size(m_NumberOfInputs,m_NumberOfInputs);
        m_InverseClassesVariances[i].set_size(m_NumberOfInputs,m_NumberOfInputs);
        m_ClassesMeans[i].resize(m_NumberOfInputs);
    }

    this->EstimateClassesParameters();

    for (unsigned int i = 0;i < m_NumberOfClasses;++i)
    {
        m_OldClassesMeans[i] = m_ClassesMeans[i];
        m_OldClassesVariances[i] = m_ClassesVariances[i];
    }

    m_LabelMap = nullptr;
    m_ZScoreMap = nullptr;
}

template <typename TInputImage>
void
TissuesEMClassificationImageFilter <TInputImage>
::GenerateData()
{
    this->AllocateOutputs();
    this->BeforeThreadedGenerateData();

    unsigned int itncount = 0;
    bool continueLoop = true;
    while ((itncount < m_MaximumIterations) && (continueLoop))
    {
        if (m_Verbose)
            std::cout << "Iteration " << itncount + 1 << "..." << std::endl;

        itk::TimeProbe tmpTime;
        tmpTime.Start();

        // Parallelize by calling DynamicThreadedGenerateData, estimation of reference
        this->GetMultiThreader()->template ParallelizeImageRegion<TOutputImage::ImageDimension> (
            this->GetComputationRegion(),
            [this](const OutputImageRegionType & outputRegionForThread)
              { this->DynamicThreadedGenerateData(outputRegionForThread); }, this);

        tmpTime.Stop();

        if (m_Verbose)
            std::cout << "Reference estimated in " << tmpTime.GetTotal() << "..." << std::endl;

        itk::TimeProbe tmpTimeParams;
        tmpTimeParams.Start();

        EstimateClassesParameters();

        tmpTimeParams.Stop();

        if (m_Verbose)
            std::cout << "Classes parameters estimated in " << tmpTimeParams.GetTotal() << "..." << std::endl;

        ++itncount;

        if (itncount != 1)
            continueLoop = !endConditionReached();

        if (continueLoop)
        {
            for (unsigned int i = 0;i < m_NumberOfClasses;++i)
            {
                m_OldClassesMeans[i] = m_ClassesMeans[i];
                m_OldClassesVariances[i] = m_ClassesVariances[i];
            }
        }
    }
}

template <typename TInputImage>
void
TissuesEMClassificationImageFilter <TInputImage>
::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    typedef itk::ImageRegionConstIterator <TInputImage> InIteratorType;
    typedef itk::ImageRegionIterator <TOutputImage> OutRegionIteratorType;
    typedef itk::ImageRegionConstIterator <LocalPriorImageType> LocalPriorIteratorType;
    typedef itk::ImageRegionIterator <MaskImageType> MaskRegionIteratorType;

    OutRegionIteratorType outItr(this->GetOutput(),outputRegionForThread);

    LocalPriorIteratorType localPriorIt(m_LocalPriorImage,outputRegionForThread);

    std::vector <InIteratorType> inputIterators(m_NumberOfInputs);
    for (unsigned int i = 0;i < m_NumberOfInputs;++i)
        inputIterators[i] = InIteratorType(this->GetInput(i), outputRegionForThread);

    MaskRegionIteratorType maskItr(this->GetComputationMask(),outputRegionForThread);

    OutputPixelType classesProbabilities(m_NumberOfClasses);
    OutputPixelType priorProbabilities(m_NumberOfClasses);
    std::vector <double> inputValues(m_NumberOfInputs);
    double piRoot = std::pow(2.0 * M_PI, m_NumberOfInputs / 2.0);

    while (!outItr.IsAtEnd())
    {
        if (maskItr.Get() == 0)
        {
            ++outItr;
            for (unsigned int i = 0;i < m_NumberOfInputs;++i)
                ++inputIterators[i];
            ++maskItr;
            ++localPriorIt;

            continue;
        }

        classesProbabilities.Fill(0.0);
        priorProbabilities = localPriorIt.Get();

        for (unsigned int i = 0;i < m_NumberOfInputs;++i)
            inputValues[i] = inputIterators[i].Get();

        double denom = 0;
        for (unsigned int m = 0;m < m_NumberOfClasses;++m)
        {
            double quadForm = 0;

            for (unsigned int i = 0;i < m_NumberOfInputs;++i)
            {
                double diff = m_ClassesMeans[m][i] - inputValues[i];
                quadForm += m_InverseClassesVariances[m](i,i) * diff * diff;
                for (unsigned int j = i+1;j < m_NumberOfInputs;++j)
                    quadForm += 2 * m_InverseClassesVariances[m](i,j) * diff * (m_ClassesMeans[m][j] - inputValues[j]);
            }

            double classValue = std::exp(- 0.5 * quadForm);
            classValue /= m_ClassesVariancesSqrtDeterminants[m] * piRoot;

            denom += priorProbabilities[m] * classValue;
            classesProbabilities[m] = classValue * priorProbabilities[m];
        }

        if (denom != 0.0)
        {
            for (unsigned int m = 0;m < m_NumberOfClasses;++m)
                classesProbabilities[m] /= denom;
        }

        outItr.Set(classesProbabilities);

        ++outItr;
        ++maskItr;
        for (unsigned int i = 0;i < m_NumberOfInputs;++i)
            ++inputIterators[i];
        ++localPriorIt;
    }
}

template <typename TInputImage>
void
TissuesEMClassificationImageFilter <TInputImage>
::EstimateClassesParameters()
{
    itk::PoolMultiThreader::Pointer threaderMstep = itk::PoolMultiThreader::New();

    EMStepThreadStruct *tmpStr = new EMStepThreadStruct;
    tmpStr->Filter = this;

    unsigned int actualNumberOfThreads = std::min(this->GetNumberOfWorkUnits(),m_NumberOfClasses);

    threaderMstep->SetNumberOfWorkUnits(actualNumberOfThreads);
    threaderMstep->SetSingleMethod(this->ThreadEstimateClassesParams,tmpStr);
    threaderMstep->SingleMethodExecute();

    delete tmpStr;
}

template <typename TInputImage>
ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION
TissuesEMClassificationImageFilter <TInputImage>
::ThreadEstimateClassesParams(void *arg)
{
    itk::MultiThreaderBase::WorkUnitInfo *threadArgs = (itk::MultiThreaderBase::WorkUnitInfo *)arg;

    unsigned int nbThread = threadArgs->WorkUnitID;
    unsigned int nbProcs = threadArgs->NumberOfWorkUnits;

    EMStepThreadStruct *tmpStr = (EMStepThreadStruct *)threadArgs->UserData;
    unsigned int numClasses = tmpStr->Filter->GetNumberOfClasses();
    unsigned int actualNumberOfThreads = std::min(nbProcs,numClasses);

    unsigned int classStart = (unsigned int)std::floor((double)nbThread*numClasses/actualNumberOfThreads);
    unsigned int classEnd = (unsigned int)std::floor((double)(nbThread + 1.0)*numClasses/actualNumberOfThreads);

    classEnd = std::min(numClasses,classEnd);

    tmpStr->Filter->EstimateClassesParameters(classStart,classEnd);

    return ITK_THREAD_RETURN_DEFAULT_VALUE;
}

template <typename TInputImage>
void
TissuesEMClassificationImageFilter <TInputImage>
::EstimateClassesParameters(unsigned int classStart, unsigned int classEnd)
{
    typedef itk::ImageRegionConstIterator <TInputImage> InIteratorType;
    typedef itk::ImageRegionIterator <TOutputImage> OutRegionIteratorType;
    typedef itk::ImageRegionIterator <MaskImageType> MaskRegionIteratorType;

    OutRegionIteratorType outItr(this->GetOutput(),this->GetComputationRegion());
    MaskRegionIteratorType maskItr(this->GetComputationMask(),this->GetComputationRegion());
    std::vector <InIteratorType> inputIterators(m_NumberOfInputs);
    for (unsigned int i = 0;i < m_NumberOfInputs;++i)
        inputIterators[i] = InIteratorType(this->GetInput(i),this->GetComputationRegion());

    OutputPixelType classesVector(m_NumberOfClasses);
    std::vector <double> inputValues(m_NumberOfInputs);
    anima::LogEuclideanTensorCalculator <double>::Pointer leCalc = anima::LogEuclideanTensorCalculator <double>::New();

    for (unsigned int i = classStart; i < classEnd;++i)
    {
        for (unsigned int j = 0;j < m_NumberOfInputs;++j)
           inputIterators[j].GoToBegin();
        outItr.GoToBegin();
        maskItr.GoToBegin();

        double classDenom = 0.0;
        for (unsigned int j = 0;j < m_NumberOfInputs;++j)
            m_ClassesMeans[i][j] = 0;
        m_ClassesVariances[i].fill(0.0);

        while (!maskItr.IsAtEnd())
        {
            if (maskItr.Get() == 0)
            {
                ++maskItr;
                ++outItr;
                for (unsigned int j = 0;j < m_NumberOfInputs;++j)
                    ++inputIterators[j];

                continue;
            }

            classesVector = outItr.Get();
            classDenom += classesVector[i];

            for (unsigned int j = 0;j < m_NumberOfInputs;++j)
            {
                double inputValue = inputIterators[j].Get();
                m_ClassesMeans[i][j] += classesVector[i] * inputValue;
            }

            ++maskItr;
            ++outItr;
            for (unsigned int j = 0;j < m_NumberOfInputs;++j)
                ++inputIterators[j];
        }

        for (unsigned int j = 0;j < m_NumberOfInputs;++j)
            m_ClassesMeans[i][j] /= classDenom;

        for (unsigned int j = 0;j < m_NumberOfInputs;++j)
            inputIterators[j].GoToBegin();
        outItr.GoToBegin();
        maskItr.GoToBegin();

        while (!maskItr.IsAtEnd())
        {
            if (maskItr.Get() == 0)
            {
                ++maskItr;
                ++outItr;
                for (unsigned int j = 0;j < m_NumberOfInputs;++j)
                    ++inputIterators[j];

                continue;
            }

            classesVector = outItr.Get();

            for (unsigned int j = 0;j < m_NumberOfInputs;++j)
                inputValues[j] = inputIterators[j].Get();

            for (unsigned int j = 0;j < m_NumberOfInputs;++j)
            {
                m_ClassesVariances[i](j,j) += classesVector[i] * (inputValues[j] - m_ClassesMeans[i][j]) * (inputValues[j] - m_ClassesMeans[i][j]);
                for (unsigned int k = j + 1;k < m_NumberOfInputs;++k)
                    m_ClassesVariances[i](j,k) += classesVector[i] * (inputValues[j] - m_ClassesMeans[i][j]) * (inputValues[k] - m_ClassesMeans[i][k]);
            }

            ++maskItr;
            ++outItr;
            for (unsigned int j = 0;j < m_NumberOfInputs;++j)
                ++inputIterators[j];
        }

        for (unsigned int j = 0;j < m_NumberOfInputs;++j)
        {
            m_ClassesVariances[i](j,j) /= classDenom;
            for (unsigned int k = j + 1;k < m_NumberOfInputs;++k)
            {
                m_ClassesVariances[i](j,k) /= classDenom;
                m_ClassesVariances[i](k,j) = m_ClassesVariances[i](j,k);
            }
        }

        m_ClassesVariancesSqrtDeterminants[i] = std::sqrt(vnl_determinant(m_ClassesVariances[i]));
        leCalc->GetTensorPower(m_ClassesVariances[i],m_InverseClassesVariances[i], -1.0);
    }
}

template <typename TInputImage>
typename TissuesEMClassificationImageFilter <TInputImage>::MaskImagePointer &
TissuesEMClassificationImageFilter <TInputImage>
::GetClassificationAsLabelMap()
{
    if (m_LabelMap)
        return m_LabelMap;

    m_LabelMap = MaskImageType::New();

    m_LabelMap->Initialize();
    m_LabelMap->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    m_LabelMap->SetSpacing(this->GetInput()->GetSpacing());
    m_LabelMap->SetOrigin(this->GetInput()->GetOrigin());
    m_LabelMap->SetDirection(this->GetInput()->GetDirection());

    m_LabelMap->Allocate();
    m_LabelMap->FillBuffer(0);

    typedef itk::ImageRegionIterator <TOutputImage> OutRegionIteratorType;
    typedef itk::ImageRegionIterator <MaskImageType> MaskRegionIteratorType;

    InputImageRegionType largestRegion = this->GetComputationRegion();
    MaskRegionIteratorType labelMapItr(m_LabelMap,largestRegion);
    OutRegionIteratorType outItr(this->GetOutput(),largestRegion);
    MaskRegionIteratorType maskItr(this->GetComputationMask(),largestRegion);

    while(!maskItr.IsAtEnd())
    {
        if (maskItr.Get() != 0)
        {
            OutputPixelType tmpVal = outItr.Get();
            unsigned int maxClass = 0;
            double maxVal = tmpVal[0];

            for (unsigned int i = 0;i < m_NumberOfClasses;++i)
            {
                if (tmpVal[i] > maxVal)
                {
                    maxVal = tmpVal[i];
                    maxClass = i;
                }
            }

            labelMapItr.Set(maxClass + 1);
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
typename TissuesEMClassificationImageFilter <TInputImage>::RealImagePointer &
TissuesEMClassificationImageFilter <TInputImage>
::GetZScoreMap()
{
    this->GetClassificationAsLabelMap();

    m_ZScoreMap = RealImageType::New();

    m_ZScoreMap->Initialize();
    m_ZScoreMap->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    m_ZScoreMap->SetSpacing(this->GetInput()->GetSpacing());
    m_ZScoreMap->SetOrigin(this->GetInput()->GetOrigin());
    m_ZScoreMap->SetDirection(this->GetInput()->GetDirection());

    m_ZScoreMap->Allocate();
    m_ZScoreMap->FillBuffer(0);

    typedef itk::ImageRegionConstIterator <TInputImage> InRegionIteratorType;
    typedef itk::ImageRegionConstIterator <MaskImageType> MaskRegionIteratorType;
    typedef itk::ImageRegionIterator <RealImageType> RealRegionIteratorType;

    InputImageRegionType largestRegion = this->GetComputationRegion();
    MaskRegionIteratorType labelMapItr(m_LabelMap,largestRegion);
    MaskRegionIteratorType maskItr(this->GetComputationMask(),largestRegion);
    RealRegionIteratorType zscItr(m_ZScoreMap,largestRegion);
    std::vector <InRegionIteratorType> inputIterators(m_NumberOfInputs);
    for (unsigned int i = 0;i < m_NumberOfInputs;++i)
        inputIterators[i] = InRegionIteratorType(this->GetInput(i),largestRegion);

    std::vector <double> inputValues(m_NumberOfInputs);
    while(!maskItr.IsAtEnd())
    {
        if (maskItr.Get() != 0)
        {
            unsigned int classSelected = labelMapItr.Get() - 1;

            double zsc = 0;

            for (unsigned int i = 0;i < m_NumberOfInputs;++i)
                inputValues[i] = inputIterators[i].Get();

            for (unsigned int i = 0;i < m_NumberOfInputs;++i)
            {
                double diff = m_ClassesMeans[classSelected][i] - inputValues[i];
                zsc += m_InverseClassesVariances[classSelected](i,i) * diff * diff;
                for (unsigned int j = i+1;j < m_NumberOfInputs;++j)
                    zsc += 2 * m_InverseClassesVariances[classSelected](i,j) * diff * (m_ClassesMeans[classSelected][j] - inputValues[j]);
            }

            zscItr.Set(std::sqrt(zsc));
        }
        else
            zscItr.Set(0);

        ++labelMapItr;
        ++maskItr;
        ++zscItr;
        for (unsigned int i = 0;i < m_NumberOfInputs;++i)
            ++inputIterators[i];
    }

    return m_ZScoreMap;
}

} // end namespace anima
