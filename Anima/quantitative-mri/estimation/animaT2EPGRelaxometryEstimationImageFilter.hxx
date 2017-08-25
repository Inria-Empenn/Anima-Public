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

    typedef anima::BobyqaOptimizer OptimizerType;
    typedef anima::T2EPGRelaxometryCostFunction CostFunctionType;

    typename CostFunctionType::Pointer cost = CostFunctionType::New();
    cost->SetT2EchoSpacing(m_EchoSpacing);
    cost->SetT2ExcitationFlipAngle(m_T2ExcitationFlipAngle);
    cost->SetT2FlipAngles(m_T2FlipAngles);
    cost->SetB1OnExcitationAngle(m_B1OnExcitationAngle);

    while (!maskItr.IsAtEnd())
    {
        double t1Value = 1;

        if (m_T1Map)
            t1Value = t1MapItr.Get();

        double b1Value = 1;
        double b1ValueOld = 0;
        double t2Value = m_T2UpperBound / 2.0;
        if (m_T1Map && (m_T2UpperBound > t1Value))
            t2Value = t1Value / 2.0;
        double m0Value = m_M0UpperBound / 2.0;

        if ((maskItr.Get() == 0)||(b1Value == 0)||(t1Value == 0))
        {
            outT2Iterator.Set(0);
            outM0Iterator.Set(0);
            outB1Iterator.Set(0);

            ++maskItr;
            ++outT2Iterator;
            ++outM0Iterator;
            ++outB1Iterator;

            for (unsigned int i = 0;i < numInputs;++i)
                ++inIterators[i];

            if (m_T1Map)
                ++t1MapItr;

            continue;
        }

        cost->SetT1Value(t1Value);
        cost->SetM0Value(m0Value);
        cost->SetT2Value(t2Value);
        cost->SetB1Value(b1Value);

        // Here go the T2 and M0 estimation
        for (unsigned int i = 0;i < numInputs;++i)
            relaxoT2Data[i] = inIterators[i].Get();

        cost->SetT2RelaxometrySignals(relaxoT2Data);

        unsigned int iterNum = 0;
        while ((std::abs(b1ValueOld - b1Value) / b1Value > m_OptimizerStopCondition)&&(iterNum < m_MaximumOptimizerIterations))
        {
            b1ValueOld = b1Value;
            ++iterNum;

            // Switch between B1 and T2 optimization
            for (unsigned int i = 0;i < 2;++i)
            {
                cost->SetOptimizeB1Value(i == 1);
                unsigned int dimension = cost->GetNumberOfParameters();
                itk::Array<double> lowerBounds(dimension);
                itk::Array<double> upperBounds(dimension);
                OptimizerType::ParametersType p(dimension);
                OptimizerType::ScalesType scales(dimension);

                OptimizerType::Pointer optimizer = OptimizerType::New();
                optimizer->SetRhoBegin(m_OptimizerInitialStep);
                optimizer->SetRhoEnd(m_OptimizerStopCondition);
                optimizer->SetNumberSamplingPoints(dimension + 2);
                optimizer->SetMaximumIteration(m_MaximumOptimizerIterations);

                if (i == 0) // T2 optimization
                {
                    lowerBounds[0] = 1.0e-4;
                    upperBounds[0] = m_T2UpperBound;
                    if (m_T1Map && (m_T2UpperBound > t1Value))
                        upperBounds[0] = t1Value;

                    lowerBounds[1] = 1.0e-4;
                    upperBounds[1] = m_M0UpperBound;

                    p[0] = t2Value;
                    p[1] = m0Value;
                    scales[0] = (upperBounds[1] - lowerBounds[1]) / (upperBounds[0] - lowerBounds[0]);
                    scales[1] = 1;
                }
                else
                {
                    lowerBounds[0] = 0.5;
                    upperBounds[0] = 2.0;
                    p[0] = b1Value;
                    scales[0] = 1.0;
                }

                optimizer->SetScales(scales);
                optimizer->SetLowerBounds(lowerBounds);
                optimizer->SetUpperBounds(upperBounds);
                optimizer->SetMaximize(false);

                optimizer->SetCostFunction(cost);

                optimizer->SetInitialPosition(p);
                optimizer->StartOptimization();

                p = optimizer->GetCurrentPosition();

                if (i == 0)
                {
                    t2Value = p[0];
                    m0Value = p[1];
                    cost->SetM0Value(m0Value);
                    cost->SetT2Value(t2Value);
                }
                else
                {
                    b1Value = p[0];
                    cost->SetB1Value(b1Value);
                }
            }
        }

        outT2Iterator.Set(t2Value);
        outM0Iterator.Set(m0Value);
        outB1Iterator.Set(b1ValueOld);

        ++maskItr;
        ++outT2Iterator;
        ++outM0Iterator;
        ++outB1Iterator;

        for (unsigned int i = 0;i < numInputs;++i)
            ++inIterators[i];

        if (m_T1Map)
            ++t1MapItr;
    }
}

} // end namespace anima
