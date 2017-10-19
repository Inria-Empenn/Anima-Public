#pragma once
#include "animaT2RelaxometryEstimationImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

#include <animaNLOPTOptimizers.h>
#include <animaT2RelaxometryCostFunction.h>

namespace anima
{

template <typename TInputImage, typename TOutputImage>
void
T2RelaxometryEstimationImageFilter <TInputImage,TOutputImage>
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
T2RelaxometryEstimationImageFilter <TInputImage,TOutputImage>
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

    typedef anima::NLOPTOptimizers OptimizerType;
    typedef anima::T2RelaxometryCostFunction CostFunctionType;

    typename CostFunctionType::Pointer cost = CostFunctionType::New();
    cost->SetT2EchoSpacing(m_EchoSpacing);
    cost->SetTRValue(m_TRValue);

    unsigned int dimension = cost->GetNumberOfParameters();
    itk::Array<double> lowerBounds(dimension);
    itk::Array<double> upperBounds(dimension);
    OptimizerType::ParametersType p(dimension);

    std::vector <double> relaxoT2Data(numInputs,0);

    while (!maskItr.IsAtEnd())
    {
        if (maskItr.Get() == 0)
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

            continue;
        }

        double t1Value = 1000.0;

        if (m_T1Map)
            t1Value = t1MapItr.Get();

        for (unsigned int i = 0;i < numInputs;++i)
            relaxoT2Data[i] = inIterators[i].Get();

        cost->SetT2RelaxometrySignals(relaxoT2Data);

        OptimizerType::Pointer optimizer = OptimizerType::New();
        optimizer->SetAlgorithm(NLOPT_LN_BOBYQA);
        optimizer->SetXTolRel(m_OptimizerStopCondition);
        optimizer->SetFTolRel(1.0e-2 * m_OptimizerStopCondition);
        optimizer->SetMaxEval(5000);
        optimizer->SetVectorStorageSize(2000);

        lowerBounds[0] = 1.0e-4;
        upperBounds[0] = m_T2UpperBoundValue;
        if (m_T1Map && (m_T2UpperBoundValue > t1Value))
            upperBounds[0] = t1Value;

        p[0] = 100;
        if (p[0] > upperBounds[0])
            p[0] = (upperBounds[0] - lowerBounds[0]) / 2.0;

        optimizer->SetLowerBoundParameters(lowerBounds);
        optimizer->SetUpperBoundParameters(upperBounds);
        optimizer->SetMaximize(false);

        optimizer->SetCostFunction(cost);

        optimizer->SetInitialPosition(p);
        optimizer->StartOptimization();

        p = optimizer->GetCurrentPosition();

        cost->GetValue(p);

        outT2Iterator.Set(p[0]);
        outM0Iterator.Set(cost->GetM0Value());

        ++maskItr;
        ++outT2Iterator;
        ++outM0Iterator;

        for (unsigned int i = 0;i < numInputs;++i)
            ++inIterators[i];

        if (m_T1Map)
            ++t1MapItr;
    }
}

} // end namespace anima
