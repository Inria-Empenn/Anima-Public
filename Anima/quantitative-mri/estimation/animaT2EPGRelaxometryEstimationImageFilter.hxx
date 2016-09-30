#pragma once
#include "animaT2EPGRelaxometryEstimationImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

#include <animaT2EPGRelaxometryCostFunction.h>
#include <animaBobyqaOptimizer.h>

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
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionConstIterator <InputImageType> ImageIteratorType;
    typedef itk::ImageRegionIterator <OutputImageType> OutImageIteratorType;

    OutImageIteratorType outT2Iterator(this->GetOutput(0),outputRegionForThread);
    OutImageIteratorType outM0Iterator(this->GetOutput(1),outputRegionForThread);

    unsigned int numInputs = this->GetNumberOfIndexedInputs();
    std::vector <ImageIteratorType> inIterators;
    for (unsigned int i = 0;i < numInputs;++i)
        inIterators.push_back(ImageIteratorType(this->GetInput(i),outputRegionForThread));

    typedef itk::ImageRegionIterator <MaskImageType> MaskIteratorType;
    MaskIteratorType maskItr(this->GetComputationMask(),outputRegionForThread);

    OutImageIteratorType t1MapItr;
    if (m_T1Map)
        t1MapItr = OutImageIteratorType(m_T1Map,outputRegionForThread);

    OutImageIteratorType b1MapItr;
    if (m_B1Map)
        b1MapItr = OutImageIteratorType(m_B1Map,outputRegionForThread);

    OutImageIteratorType initialT2Itr;
    if (m_InitialT2Map)
        initialT2Itr = OutImageIteratorType(m_InitialT2Map,outputRegionForThread);

    OutImageIteratorType initialM0Itr;
    if (m_InitialM0Map)
        initialM0Itr = OutImageIteratorType(m_InitialM0Map,outputRegionForThread);

    std::vector <double> relaxoT2Data(numInputs,0);

    typedef anima::BobyqaOptimizer OptimizerType;
    typedef anima::T2EPGRelaxometryCostFunction CostFunctionType;

    typename CostFunctionType::Pointer cost = CostFunctionType::New();
    cost->SetT2EchoSpacing(m_EchoSpacing);
    cost->SetT2ExcitationFlipAngle(m_T2ExcitationFlipAngle);
    cost->SetT2FlipAngles(m_T2FlipAngles);
    cost->SetB1OnExcitationAngle(m_B1OnExcitationAngle);

    unsigned int dimension = cost->GetNumberOfParameters();
    itk::Array<double> lowerBounds(dimension);
    itk::Array<double> upperBounds(dimension);
    OptimizerType::ParametersType p(dimension);
    OptimizerType::ScalesType scales(dimension);

    while (!maskItr.IsAtEnd())
    {
        double t1Value = 1;

        if (m_T1Map)
            t1Value = t1MapItr.Get();

        double b1Value = 1;

        if (m_B1Map)
            b1Value = b1MapItr.Get();

        if ((maskItr.Get() == 0)||(b1Value == 0)||(t1Value == 0))
        {
            outT2Iterator.Set(0);
            outM0Iterator.Set(0);

            ++maskItr;
            ++outT2Iterator;
            ++outM0Iterator;

            for (unsigned int i = 0;i < numInputs;++i)
                ++inIterators[i];

            if (m_T1Map)
                ++t1MapItr;

            if (m_B1Map)
                ++b1MapItr;

            if (m_InitialT2Map)
                ++initialT2Itr;

            if (m_InitialM0Map)
                ++initialM0Itr;

            continue;
        }

        cost->SetT1Value(t1Value);
        cost->SetB1Value(b1Value);

        // Here go the T2 and M0 estimation
        for (unsigned int i = 0;i < numInputs;++i)
            relaxoT2Data[i] = inIterators[i].Get();

        cost->SetT2RelaxometrySignals(relaxoT2Data);

        OptimizerType::Pointer optimizer = OptimizerType::New();
        optimizer->SetRhoBegin(m_OptimizerInitialStep);
        optimizer->SetRhoEnd(m_OptimizerStopCondition);
        optimizer->SetNumberSamplingPoints(dimension + 2);
        optimizer->SetMaximumIteration(m_MaximumOptimizerIterations);

        lowerBounds[0] = 1.0e-4;
        upperBounds[0] = m_T2UpperBound;
        if (m_T1Map && (m_T2UpperBound > t1Value))
            upperBounds[0] = t1Value;

        lowerBounds[1] = 1.0e-4;
        upperBounds[1] = m_M0UpperBound;

        if (m_InitialT2Map)
            p[0] = initialT2Itr.Get();
        else
            p[0] = (lowerBounds[0] + upperBounds[0]) / 2.0;

        if (m_InitialM0Map)
            p[1] = initialM0Itr.Get();
        else
            p[1] = (lowerBounds[1] + upperBounds[1]) / 2.0;

        scales[0] = (upperBounds[1] - lowerBounds[1]) / (upperBounds[0] - lowerBounds[0]);
        scales[1] = 1;

        optimizer->SetScales(scales);
        optimizer->SetLowerBounds(lowerBounds);
        optimizer->SetUpperBounds(upperBounds);
        optimizer->SetMaximize(false);

        optimizer->SetCostFunction(cost);

        optimizer->SetInitialPosition(p);
        optimizer->StartOptimization();

        p = optimizer->GetCurrentPosition();

        outT2Iterator.Set(p[0]);
        outM0Iterator.Set(p[1]);

        ++maskItr;
        ++outT2Iterator;
        ++outM0Iterator;

        for (unsigned int i = 0;i < numInputs;++i)
            ++inIterators[i];

        if (m_T1Map)
            ++t1MapItr;

        if (m_B1Map)
            ++b1MapItr;

        if (m_InitialT2Map)
            ++initialT2Itr;

        if (m_InitialM0Map)
            ++initialM0Itr;
    }
}

} // end namespace anima
