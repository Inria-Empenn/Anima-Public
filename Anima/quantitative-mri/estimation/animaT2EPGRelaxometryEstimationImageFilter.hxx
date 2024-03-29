#pragma once
#include "animaT2EPGRelaxometryEstimationImageFilter.h"
#include <animaT2RelaxometryEstimationImageFilter.h>

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

#include <animaT2EPGRelaxometryCostFunction.h>
#include <animaNLOPTOptimizers.h>

namespace anima
{

template <typename TInputImage, typename TOutputImage>
void
T2EPGRelaxometryEstimationImageFilter <TInputImage,TOutputImage>
::CheckComputationMask()
{
    if (this->GetComputationMask())
        return;

    typedef itk::ImageRegionConstIterator <InputImageType> IteratorType;
    typedef itk::ImageRegionIterator <MaskImageType> MaskIteratorType;

    std::vector <IteratorType> inItrs(this->GetNumberOfIndexedInputs());
    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
        inItrs[i] = IteratorType(this->GetInput(i),this->GetOutput()->GetLargestPossibleRegion());

    typename MaskImageType::Pointer maskImage = MaskImageType::New();
    maskImage->Initialize();
    maskImage->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
    maskImage->SetSpacing (this->GetInput(0)->GetSpacing());
    maskImage->SetOrigin (this->GetInput(0)->GetOrigin());
    maskImage->SetDirection (this->GetInput(0)->GetDirection());
    maskImage->Allocate();

    MaskIteratorType maskItr (maskImage,this->GetOutput()->GetLargestPossibleRegion());
    while (!maskItr.IsAtEnd())
    {
        double averageVal = 0;
        for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
            averageVal += inItrs[i].Get();

        averageVal /= this->GetNumberOfIndexedInputs();

        bool maskPoint = (averageVal <= m_AverageSignalThreshold);

        if (maskPoint)
            maskItr.Set(0);
        else
            maskItr.Set(1);

        for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
            ++inItrs[i];

        ++maskItr;
    }

    this->SetComputationMask(maskImage);
}

template <typename TInputImage, typename TOutputImage>
void
T2EPGRelaxometryEstimationImageFilter <TInputImage,TOutputImage>
::BeforeThreadedGenerateData()
{
    this->Superclass::BeforeThreadedGenerateData();

    typedef anima::T2RelaxometryEstimationImageFilter <TInputImage,TOutputImage> InitT2FilterType;
    typename InitT2FilterType::Pointer initFilter = InitT2FilterType::New();

    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
        initFilter->SetInput(i,this->GetInput(i));

    initFilter->SetT1Map(m_T1Map);

    initFilter->SetEchoSpacing(m_EchoSpacing);
    initFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    initFilter->SetComputationMask(this->GetComputationMask());
    initFilter->SetTRValue(m_TRValue);
    initFilter->SetT2UpperBoundValue(m_T2UpperBound);
    initFilter->SetOptimizerStopCondition(m_OptimizerStopCondition);
    initFilter->SetMaximumOptimizerIterations(m_MaximumOptimizerIterations);

    initFilter->Update();

    m_InitialT2Image = initFilter->GetOutput();
    m_InitialT2Image->DisconnectPipeline();

    // Prepare pixel width data from pulse profiles if present
    if (!m_UniformPulses)
    {
        double sliceThickness = this->GetInput(0)->GetSpacing()[2];
        double excitationRatioToProfile = sliceThickness / m_ReferenceSliceThickness;
        for (unsigned int i = 0;i < m_ExcitationProfile.size();++i)
            m_ExcitationProfile[i].first *= excitationRatioToProfile;

        unsigned int zeroIndex = 0;
        while (m_ExcitationProfile[zeroIndex].first < 0)
            ++zeroIndex;

        bool continueLoop = true;
        while (continueLoop)
        {
            ++zeroIndex;
            if (zeroIndex == m_ExcitationProfile.size())
            {
                continueLoop = false;
                break;
            }

            double diff = m_ExcitationProfile[zeroIndex].second - m_ExcitationProfile[zeroIndex - 1].second;
            if ((diff > 0) && (m_ExcitationProfile[zeroIndex - 1].second < 0.01))
            {
                --zeroIndex;
                continueLoop = false;
            }
        }

        m_ExcitationPixelWidth = m_ExcitationProfile[zeroIndex].first * 2.0;

        double pulseRatioToProfile = sliceThickness * m_PulseWidthFactor / m_ReferenceSliceThickness;
        for (unsigned int i = 0;i < m_PulseProfile.size();++i)
            m_PulseProfile[i].first *= pulseRatioToProfile;
    }
}

template <typename TInputImage, typename TOutputImage>
void
T2EPGRelaxometryEstimationImageFilter <TInputImage,TOutputImage>
::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    typedef itk::ImageRegionConstIterator <InputImageType> ImageIteratorType;
    typedef itk::ImageRegionIterator <OutputImageType> OutImageIteratorType;

    OutImageIteratorType initT2Iterator(m_InitialT2Image,outputRegionForThread);
    OutImageIteratorType outT2Iterator(this->GetOutput(0),outputRegionForThread);
    OutImageIteratorType outM0Iterator(this->GetOutput(1),outputRegionForThread);
    OutImageIteratorType outB1Iterator(this->GetOutput(2),outputRegionForThread);

    unsigned int numInputs = this->GetNumberOfIndexedInputs();
    std::vector <ImageIteratorType> inIterators;
    for (unsigned int i = 0;i < numInputs;++i)
        inIterators.push_back(ImageIteratorType(this->GetInput(i),outputRegionForThread));

    typedef itk::ImageRegionIterator <MaskImageType> MaskIteratorType;
    MaskIteratorType maskItr(this->GetComputationMask(),outputRegionForThread);

    OutImageIteratorType t1MapItr;
    if (m_T1Map)
        t1MapItr = OutImageIteratorType(m_T1Map,outputRegionForThread);

    std::vector <double> relaxoT2Data(numInputs,0);

    typedef anima::NLOPTOptimizers OptimizerType;
    typedef anima::T2EPGRelaxometryCostFunction CostFunctionType;

    typename CostFunctionType::Pointer cost = CostFunctionType::New();
    cost->SetT2EchoSpacing(m_EchoSpacing);
    cost->SetT2ExcitationFlipAngle(m_T2ExcitationFlipAngle);
    cost->SetT2FlipAngles(m_T2FlipAngles);

    cost->SetUniformPulses(m_UniformPulses);
    if (!m_UniformPulses)
    {
        cost->SetPulseProfile(m_PulseProfile);
        cost->SetExcitationProfile(m_ExcitationProfile);
        cost->SetPixelWidth(m_ExcitationPixelWidth);
    }

    unsigned int dimension = cost->GetNumberOfParameters();
    itk::Array<double> lowerBounds(dimension);
    itk::Array<double> upperBounds(dimension);
    OptimizerType::ParametersType p(dimension);

    while (!maskItr.IsAtEnd())
    {
        double t1Value = m_T2UpperBound;

        if (m_T1Map)
        {
            t1Value = t1MapItr.Get();
            if (t1Value <= 0.0)
                t1Value = 1000;
        }

        double b1Value = 0.9;
        double t2Value = initT2Iterator.Get();

        if (maskItr.Get() == 0)
        {
            outT2Iterator.Set(0);
            outM0Iterator.Set(0);
            outB1Iterator.Set(0);

            ++maskItr;
            ++outT2Iterator;
            ++outM0Iterator;
            ++outB1Iterator;
            ++initT2Iterator;

            for (unsigned int i = 0;i < numInputs;++i)
                ++inIterators[i];

            if (m_T1Map)
                ++t1MapItr;

            continue;
        }

        cost->SetT1Value(t1Value);

        // Here go the T2 and B1 estimation
        for (unsigned int i = 0;i < numInputs;++i)
            relaxoT2Data[i] = inIterators[i].Get();

        cost->SetT2RelaxometrySignals(relaxoT2Data);

        OptimizerType::Pointer optimizer = OptimizerType::New();
        optimizer->SetAlgorithm(NLOPT_LN_BOBYQA);
        optimizer->SetXTolRel(m_OptimizerStopCondition);
        optimizer->SetFTolRel(1.0e-2 * m_OptimizerStopCondition);
        optimizer->SetMaxEval(m_MaximumOptimizerIterations);
        optimizer->SetVectorStorageSize(2000);

        lowerBounds[0] = 1.0;
        upperBounds[0] = m_T2UpperBound;
        if (m_T1Map && (m_T2UpperBound > t1Value))
            upperBounds[0] = t1Value;

        lowerBounds[1] = 0.5;
        upperBounds[1] = 1.0;

        p[0] = t2Value;
        p[1] = b1Value;

        optimizer->SetLowerBoundParameters(lowerBounds);
        optimizer->SetUpperBoundParameters(upperBounds);
        optimizer->SetMaximize(false);

        optimizer->SetCostFunction(cost);

        optimizer->SetInitialPosition(p);
        optimizer->StartOptimization();

        p = optimizer->GetCurrentPosition();

        t2Value = p[0];
        b1Value = p[1];
        cost->GetValue(p);

        outT2Iterator.Set(t2Value);
        outM0Iterator.Set(cost->GetM0Value());
        outB1Iterator.Set(b1Value);

        this->IncrementNumberOfProcessedPoints();
        ++maskItr;
        ++outT2Iterator;
        ++outM0Iterator;
        ++outB1Iterator;
        ++initT2Iterator;

        for (unsigned int i = 0;i < numInputs;++i)
            ++inIterators[i];

        if (m_T1Map)
            ++t1MapItr;
    }
}

} // end namespace anima
