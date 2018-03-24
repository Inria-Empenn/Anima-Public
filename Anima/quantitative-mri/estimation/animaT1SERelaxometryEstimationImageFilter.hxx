#pragma once
#include "animaT1SERelaxometryEstimationImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

#include <animaNLOPTOptimizers.h>
#include <animaT1SERelaxometryCostFunction.h>

namespace anima
{
    
template <typename TInputImage, typename TOutputImage>
void
T1SERelaxometryEstimationImageFilter <TInputImage,TOutputImage>
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
T1SERelaxometryEstimationImageFilter <TInputImage,TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionConstIterator <InputImageType> ImageIteratorType;
    typedef itk::ImageRegionIterator <OutputImageType> OutImageIteratorType;

    OutImageIteratorType outT1Iterator(this->GetOutput(1),outputRegionForThread);
    OutImageIteratorType outM0Iterator(this->GetOutput(0),outputRegionForThread);

    unsigned int numInputs = this->GetNumberOfIndexedInputs();
    std::vector <ImageIteratorType> inIterators;
    for (unsigned int i = 0;i < numInputs;++i)
        inIterators.push_back(ImageIteratorType(this->GetInput(i),outputRegionForThread));

    typedef itk::ImageRegionIterator <MaskImageType> MaskIteratorType;
    MaskIteratorType maskItr(this->GetComputationMask(),outputRegionForThread);

    anima::T1SERelaxometryCostFunction::Pointer cost = anima::T1SERelaxometryCostFunction::New();
    cost->SetTRValues(m_TRValues);

    std::vector <double> relaxoData(numInputs,0);
    typedef anima::NLOPTOptimizers OptimizerType;

    unsigned int dimension = cost->GetNumberOfParameters();
    itk::Array<double> lowerBounds(dimension);
    itk::Array<double> upperBounds(dimension);
    OptimizerType::ParametersType p(dimension);

    lowerBounds[0] = 1.0e-4;
    upperBounds[0] = m_M0UpperBound;

    lowerBounds[1] = 1.0e-4;
    upperBounds[1] = m_T1UpperBound;

    while (!maskItr.IsAtEnd())
    {
        if (maskItr.Get() == 0)
        {
            outT1Iterator.Set(0);
            outM0Iterator.Set(0);

            ++maskItr;
            ++outT1Iterator;
            ++outM0Iterator;

            for (unsigned int i = 0;i < numInputs;++i)
                ++inIterators[i];

            continue;
        }

        for (unsigned int i = 0;i < numInputs;++i)
            relaxoData[i] = inIterators[i].Get();

        cost->SetRelaxometrySignals(relaxoData);

        OptimizerType::Pointer optimizer = OptimizerType::New();
        optimizer->SetAlgorithm(NLOPT_LN_BOBYQA);

        optimizer->SetMaxEval(m_MaximumOptimizerIterations);
        optimizer->SetXTolRel(m_OptimizerStopCondition);

        optimizer->SetLowerBoundParameters(lowerBounds);
        optimizer->SetUpperBoundParameters(upperBounds);

        p[0] = 1500;
        p[1] = 1500;

        optimizer->SetMaximize(false);

        optimizer->SetCostFunction(cost);

        optimizer->SetInitialPosition(p);
        optimizer->StartOptimization();

        p = optimizer->GetCurrentPosition();

        outM0Iterator.Set(p[0]);
        outT1Iterator.Set(p[1]);

        ++maskItr;
        ++outT1Iterator;
        ++outM0Iterator;

        for (unsigned int i = 0;i < numInputs;++i)
            ++inIterators[i];
    }
}

} // end of namespace anima
