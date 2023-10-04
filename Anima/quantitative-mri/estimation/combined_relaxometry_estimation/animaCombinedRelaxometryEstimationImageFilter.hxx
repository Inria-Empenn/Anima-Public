#pragma once
#include "animaCombinedRelaxometryEstimationImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

#include <animaCombinedRelaxometryCostFunction.h>
#include <animaNLOPTOptimizers.h>

#include <animaT1RelaxometryEstimationImageFilter.h>
#include <animaT2RelaxometryEstimationImageFilter.h>

#include <itkMultiplyImageFilter.h>
#include <animaSmoothingRecursiveYvvGaussianImageFilter.h>

#include <itkTimeProbe.h>

namespace anima
{

template <typename TInputImage, typename TOutputImage>
void
CombinedRelaxometryEstimationImageFilter <TInputImage,TOutputImage>
::AddT1RelaxometryInput(TInputImage *image)
{
    m_T1RelaxometryInputs.push_back(image);
}

template <typename TInputImage, typename TOutputImage>
void
CombinedRelaxometryEstimationImageFilter <TInputImage,TOutputImage>
::CheckComputationMask()
{
    if (this->GetComputationMask())
        return;

    typedef itk::ImageRegionConstIterator <InputImageType> IteratorType;
    typedef itk::ImageRegionIterator <MaskImageType> MaskIteratorType;

    std::vector <IteratorType> inT2Itrs(this->GetNumberOfIndexedInputs());
    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
        inT2Itrs[i] = IteratorType(this->GetInput(i),this->GetOutput()->GetLargestPossibleRegion());

    std::vector <IteratorType> inT1Itrs(m_T1RelaxometryInputs.size());
    for (unsigned int i = 0;i < m_T1RelaxometryInputs.size();++i)
        inT1Itrs[i] = IteratorType(m_T1RelaxometryInputs[i],this->GetOutput()->GetLargestPossibleRegion());

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
            averageVal += inT2Itrs[i].Get();

        averageVal /= this->GetNumberOfIndexedInputs();

        bool maskPoint = (averageVal <= m_AverageSignalThreshold);

        averageVal = 0;
        for (unsigned int i = 0;i < m_T1RelaxometryInputs.size();++i)
            averageVal += inT1Itrs[i].Get();

        averageVal /= m_T1RelaxometryInputs.size();
        maskPoint = maskPoint || (averageVal <= m_AverageSignalThreshold);

        if (maskPoint)
            maskItr.Set(0);
        else
            maskItr.Set(1);

        for (unsigned int i = 0;i < m_T1RelaxometryInputs.size();++i)
            ++inT1Itrs[i];

        for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
            ++inT2Itrs[i];

        ++maskItr;
    }

    this->SetComputationMask(maskImage);
}

template <typename TInputImage, typename TOutputImage>
void
CombinedRelaxometryEstimationImageFilter <TInputImage,TOutputImage>
::GenerateData()
{
    this->AllocateOutputs();
    this->BeforeThreadedGenerateData();

    // Here comes the main part now
    this->GetOutput(0)->FillBuffer(0);
    this->GetOutput(1)->FillBuffer(0);
    this->GetOutput(2)->FillBuffer(0);
    this->GetOutput(3)->FillBuffer(1);
    this->GetOutput(4)->FillBuffer(0);

    // Compute initial values for all outputs
    // Threaded initialization
    this->InitializeOutputs();

    // Loop starts here
    for (unsigned int i = 0;i < m_NumberOfIterations;++i)
    {
        itk::TimeProbe tmpTimer;
        tmpTimer.Start();

        std::cout << "Performing iteration " << i << "... " << std::flush;

        // Estimate k factor for M0
        this->UpdateT1KFactorM0();

        // Then call alternate optimization threaded part
        this->GetMultiThreader()->template ParallelizeImageRegion<TOutputImage::ImageDimension> (
            this->GetOutput()->GetRequestedRegion(),
            [this](const OutputImageRegionType & outputRegionForThread)
              { this->DynamicThreadedGenerateData(outputRegionForThread); }, this);

        // Then smooth B1 map
        if (m_B1SmoothingSigma > 0)
            this->SmoothB1Field();

        tmpTimer.Stop();

        std::cout << "Done in " << tmpTimer.GetTotal() << "s" << std::endl;
    }
}

template <typename TInputImage, typename TOutputImage>
void
CombinedRelaxometryEstimationImageFilter <TInputImage,TOutputImage>
::InitializeOutputs()
{
    typedef anima::T1RelaxometryEstimationImageFilter <InputImageType,OutputImageType> T1EstimatorType;

    // T1 estimation
    typename T1EstimatorType::Pointer t1Estimator = T1EstimatorType::New();
    t1Estimator->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    t1Estimator->SetInput(0,m_T1RelaxometryInputs[0]);
    t1Estimator->SetInput(1,m_T1RelaxometryInputs[m_T1RelaxometryInputs.size()-1]);
    t1Estimator->SetTRValue(m_TRT1Value);
    t1Estimator->SetT1UpperBoundValue(m_T1UpperBound);
    t1Estimator->SetM0UpperBoundValue(m_M0UpperBound);
    t1Estimator->SetComputationMask(this->GetComputationMask());
    t1Estimator->SetFlipAngles(m_T1FlipAngles);

    t1Estimator->Update();

    OutputImagePointer tmpOutput = t1Estimator->GetOutput(0);
    tmpOutput->DisconnectPipeline();

    this->SetNthOutput(0,tmpOutput);

    // T2 estimation
    typedef anima::T2RelaxometryEstimationImageFilter <InputImageType,OutputImageType> T2EstimatorType;

    typename T2EstimatorType::Pointer t2Estimator = T2EstimatorType::New();
    t2Estimator->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
        t2Estimator->SetInput(i,this->GetInput(i));

    t2Estimator->SetTRValue(m_TRT2Value);
    t2Estimator->SetT2UpperBoundValue(m_T2UpperBound);
    t2Estimator->SetComputationMask(this->GetComputationMask());

    t2Estimator->SetEchoSpacing(m_T2EchoSpacing);
    t2Estimator->SetT1Map(this->GetOutput(0));

    t2Estimator->Update();

    tmpOutput = t2Estimator->GetOutput(0);
    tmpOutput->DisconnectPipeline();
    this->SetNthOutput(1,tmpOutput);

    tmpOutput = t2Estimator->GetOutput(1);
    tmpOutput->DisconnectPipeline();
    this->SetNthOutput(2,tmpOutput);

    // Now set values of B1 map that are in mask to 1, in case they are no;
    typedef itk::ImageRegionIterator <OutputImageType> OutImageIteratorType;
    OutImageIteratorType b1Itr(this->GetOutput(3),this->GetOutput(0)->GetLargestPossibleRegion());

    typedef itk::ImageRegionIterator <MaskImageType> MaskIteratorType;
    MaskIteratorType maskItr(this->GetComputationMask(),this->GetOutput(0)->GetLargestPossibleRegion());

    while (!maskItr.IsAtEnd())
    {
        if (maskItr.Get() != 0)
            b1Itr.Set(1);

        ++b1Itr;
        ++maskItr;
    }
}

template <typename TInputImage, typename TOutputImage>
void
CombinedRelaxometryEstimationImageFilter <TInputImage,TOutputImage>
::UpdateT1KFactorM0()
{
    typedef itk::ImageRegionConstIterator <InputImageType> ImageIteratorType;

    unsigned int numT1Inputs = m_T1RelaxometryInputs.size();
    std::vector <ImageIteratorType> inT1Iterators;
    for (unsigned int i = 0;i < numT1Inputs;++i)
        inT1Iterators.push_back(ImageIteratorType(m_T1RelaxometryInputs[i],this->GetOutput(0)->GetLargestPossibleRegion()));

    typedef itk::ImageRegionIterator <OutputImageType> OutImageIteratorType;
    OutImageIteratorType outT1Iterator(this->GetOutput(0),this->GetOutput(0)->GetLargestPossibleRegion());
    OutImageIteratorType outM0Iterator(this->GetOutput(2),this->GetOutput(0)->GetLargestPossibleRegion());
    OutImageIteratorType outB1Iterator(this->GetOutput(3),this->GetOutput(0)->GetLargestPossibleRegion());

    typedef itk::ImageRegionIterator <MaskImageType> MaskIteratorType;
    MaskIteratorType maskItr(this->GetComputationMask(),this->GetOutput(0)->GetLargestPossibleRegion());

    std::vector <double> trueValues, simulatedValues;

    while (!maskItr.IsAtEnd())
    {
        if (maskItr.Get() == 0)
        {
            for (unsigned int i = 0;i < numT1Inputs;++i)
                ++inT1Iterators[i];

            ++maskItr;
            ++outT1Iterator;
            ++outM0Iterator;
            ++outB1Iterator;

            continue;
        }

        double b1Value = outB1Iterator.Get();
        double t1Value = outT1Iterator.Get();
        double m0Value = outM0Iterator.Get();

        for (unsigned int i = 0;i < numT1Inputs;++i)
        {
            double simulatedValue = m0Value * (1.0 - std::exp(- m_TRT1Value / t1Value)) * std::sin(b1Value * m_T1FlipAngles[i]);
            simulatedValue /= 1.0 - std::exp(- m_TRT1Value / t1Value) * std::cos(b1Value * m_T1FlipAngles[i]);

            simulatedValues.push_back(simulatedValue);
            trueValues.push_back(inT1Iterators[i].Get());
         }

        for (unsigned int i = 0;i < numT1Inputs;++i)
            ++inT1Iterators[i];

        ++maskItr;
        ++outT1Iterator;
        ++outM0Iterator;
        ++outB1Iterator;
    }

    unsigned int numValues = simulatedValues.size();

    unsigned int numIter = 0;
    unsigned int numMaxIter = 100;

    std::vector <double> residuals(numValues,0.0);
    std::vector <double> weights(numValues,1.0);
    std::vector <double> leverageValues(numValues,0.0);

    double newKValue = 1.0;
    double oldKValue = 0;

    // Compute leverage values for robust fit (cf matlab robustfit)
    double meanXValue = 0;
    for (unsigned int i = 0;i < numValues;++i)
        meanXValue += simulatedValues[i];
    meanXValue /= numValues;

    double sumDiffSquared = 0;
    for (unsigned int i = 0;i < numValues;++i)
    {
        leverageValues[i] = (simulatedValues[i] - meanXValue) * (simulatedValues[i] - meanXValue);
        sumDiffSquared += leverageValues[i];
    }

    for (unsigned int i = 0;i < numValues;++i)
        leverageValues[i] /= sumDiffSquared;

    while (numIter < numMaxIter)
    {
        ++numIter;

        oldKValue = newKValue;

        double numeratorValue = 0;
        double denominatorValue = 0;

        for (unsigned int i = 0;i < numValues;++i)
        {
            numeratorValue += weights[i] * trueValues[i] * simulatedValues[i];
            denominatorValue += weights[i] * simulatedValues[i] * simulatedValues[i];
        }

        newKValue = numeratorValue / denominatorValue;

        if ((std::abs(newKValue - oldKValue) < 1.0e-4)||(!m_KFactorMEstimation))
            break;

        for (unsigned int i = 0;i < numValues;++i)
            residuals[i] = trueValues[i] - newKValue * simulatedValues[i];

        std::vector <double> sortedResiduals = residuals;
        unsigned int posMedian = floor(numValues / 2.0);
        std::partial_sort(sortedResiduals.begin(),sortedResiduals.begin() + posMedian + 1,sortedResiduals.end());

        double medianResidual = sortedResiduals[posMedian];
        for (unsigned int i = 0;i < numValues;++i)
            sortedResiduals[i] = std::abs(sortedResiduals[i] - medianResidual);
        std::partial_sort(sortedResiduals.begin(),sortedResiduals.begin() + posMedian + 1,sortedResiduals.end());

        double robustDeviation = sortedResiduals[posMedian] / 0.6745;

        for (unsigned int i = 0;i < numValues;++i)
            residuals[i] /= 4.685 * robustDeviation * std::sqrt(1 - leverageValues[i]);

        // Tukey bi-square as in matlab robust-fit
        for (unsigned int i = 0;i < numValues;++i)
        {
            weights[i] = 0;

            if (std::abs(residuals[i]) < 1.0)
                weights[i] = (1 - residuals[i] * residuals[i]) * (1 - residuals[i] * residuals[i]);
        }
    }

    m_KFactorM0 = newKValue;
}

template <typename TInputImage, typename TOutputImage>
void
CombinedRelaxometryEstimationImageFilter <TInputImage,TOutputImage>
::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    typedef itk::ImageRegionConstIterator <InputImageType> ImageIteratorType;

    unsigned int numT1Inputs = m_T1RelaxometryInputs.size();
    std::vector <ImageIteratorType> inT1Iterators;
    for (unsigned int i = 0;i < numT1Inputs;++i)
        inT1Iterators.push_back(ImageIteratorType(m_T1RelaxometryInputs[i],outputRegionForThread));

    unsigned int numT2Inputs = this->GetNumberOfIndexedInputs();
    std::vector <ImageIteratorType> inT2Iterators;
    for (unsigned int i = 0;i < numT2Inputs;++i)
        inT2Iterators.push_back(ImageIteratorType(this->GetInput(i),outputRegionForThread));

    typedef itk::ImageRegionIterator <OutputImageType> OutImageIteratorType;
    OutImageIteratorType outT1Iterator(this->GetOutput(0),outputRegionForThread);
    OutImageIteratorType outT2Iterator(this->GetOutput(1),outputRegionForThread);
    OutImageIteratorType outM0Iterator(this->GetOutput(2),outputRegionForThread);
    OutImageIteratorType outB1Iterator(this->GetOutput(3),outputRegionForThread);
    OutImageIteratorType outB1AdditiveIterator(this->GetOutput(4),outputRegionForThread);

    typedef itk::ImageRegionIterator <MaskImageType> MaskIteratorType;
    MaskIteratorType maskItr(this->GetComputationMask(),outputRegionForThread);

    std::vector <double> relaxoT1Data(numT1Inputs,0);
    std::vector <double> relaxoT2Data(numT2Inputs,0);

    typedef anima::NLOPTOptimizers OptimizerType;
    typedef anima::CombinedRelaxometryCostFunction CombinedCostFunctionType;

    CombinedCostFunctionType::Pointer cost = CombinedCostFunctionType::New();
    cost->SetT2EchoSpacing(m_T2EchoSpacing);
    cost->SetT2ExcitationFlipAngle(m_T2ExcitationFlipAngle);
    cost->SetTRValue(m_TRT1Value);
    cost->SetKValue(m_KFactorM0);
    cost->SetT1FlipAngles(m_T1FlipAngles);
    cost->SetT2FlipAngles(m_T2FlipAngles);

    unsigned int dimension = 1;
    itk::Array<double> lowerBounds(dimension);
    itk::Array<double> upperBounds(dimension);
    OptimizerType::ParametersType p(dimension);

    while (!maskItr.IsAtEnd())
    {
        if (maskItr.Get() == 0)
        {
            outT1Iterator.Set(0.0);
            outT2Iterator.Set(0.0);
            outM0Iterator.Set(0.0);

            for (unsigned int i = 0;i < numT1Inputs;++i)
                ++inT1Iterators[i];

            for (unsigned int i = 0;i < numT2Inputs;++i)
                ++inT2Iterators[i];

            ++maskItr;
            ++outT1Iterator;
            ++outT2Iterator;
            ++outM0Iterator;
            ++outB1Iterator;
            ++outB1AdditiveIterator;

            continue;
        }

        for (unsigned int i = 0;i < numT1Inputs;++i)
            relaxoT1Data[i] = inT1Iterators[i].Get();

        for (unsigned int i = 0;i < numT2Inputs;++i)
            relaxoT2Data[i] = inT2Iterators[i].Get();

        cost->SetT1RelaxometrySignals(relaxoT1Data);
        cost->SetT2RelaxometrySignals(relaxoT2Data);

        // Switch over optimizations: B1, B1 additive, T2, M0, T1
        for (unsigned int i = 0;i < 5;++i)
        {
            cost->SetT1Value(outT1Iterator.Get());
            cost->SetT2Value(outT2Iterator.Get());
            cost->SetM0Value(outM0Iterator.Get());
            cost->SetB1Value(outB1Iterator.Get());
            cost->SetB1T2AdditiveValue(outB1AdditiveIterator.Get());

            CombinedCostFunctionType::OptimizedValueType optValueType = (CombinedCostFunctionType::OptimizedValueType)i;
            cost->SetOptimizedValue(optValueType);

            OptimizerType::Pointer optimizer = OptimizerType::New();
            optimizer->SetAlgorithm(NLOPT_LN_BOBYQA);
            optimizer->SetXTolRel(m_OptimizerStopCondition);
            optimizer->SetFTolRel(1.0e-2 * m_OptimizerStopCondition);
            optimizer->SetMaxEval(m_MaximumOptimizerIterations);
            optimizer->SetVectorStorageSize(2000);

            switch(optValueType)
            {
                case CombinedCostFunctionType::B1:
                    lowerBounds[0] = m_B1LowerBound;
                    upperBounds[0] = m_B1UpperBound;
                    p[0] = outB1Iterator.Get();
                    break;

                case CombinedCostFunctionType::B1_Additive:
                {
                    double b1Value = outB1Iterator.Get();
                    p[0] = outB1AdditiveIterator.Get();

                    lowerBounds[0] = m_B1LowerBound - b1Value;
                    upperBounds[0] = 1.0 - b1Value;
                    break;
                }

                case CombinedCostFunctionType::T2:
                    lowerBounds[0] = 1.0e-4;
                    upperBounds[0] = outT1Iterator.Get();
                    if (upperBounds[0] > m_T2UpperBound)
                        upperBounds[0] = m_T2UpperBound;
                    p[0] = outT2Iterator.Get();
                    break;

                case CombinedCostFunctionType::M0:
                    lowerBounds[0] = 1.0e-4;
                    upperBounds[0] = m_M0UpperBound;
                    p[0] = outM0Iterator.Get();
                    break;

                case CombinedCostFunctionType::T1:
                default:
                    lowerBounds[0] = 1.0e-4;
                    upperBounds[0] = m_T1UpperBound;
                    p[0] = outT1Iterator.Get();
                    break;
            }

            optimizer->SetLowerBoundParameters(lowerBounds);
            optimizer->SetUpperBoundParameters(upperBounds);
            optimizer->SetMaximize(false);

            optimizer->SetCostFunction(cost);

            optimizer->SetInitialPosition(p);
            optimizer->StartOptimization();

            p = optimizer->GetCurrentPosition();

            switch(optValueType)
            {
                case CombinedCostFunctionType::B1:
                    outB1Iterator.Set(p[0]);
                    break;

                case CombinedCostFunctionType::B1_Additive:
                    outB1AdditiveIterator.Set(p[0]);
                    break;

                case CombinedCostFunctionType::T2:
                    outT2Iterator.Set(p[0]);
                    break;

                case CombinedCostFunctionType::M0:
                    outM0Iterator.Set(p[0]);
                    break;

                case CombinedCostFunctionType::T1:
                default:
                    outT1Iterator.Set(p[0]);
                    break;
            }
        }

        // Update iterators for next pixel
        for (unsigned int i = 0;i < numT1Inputs;++i)
            ++inT1Iterators[i];

        for (unsigned int i = 0;i < numT2Inputs;++i)
            ++inT2Iterators[i];

        ++maskItr;
        ++outT1Iterator;
        ++outT2Iterator;
        ++outM0Iterator;
        ++outB1Iterator;
        ++outB1AdditiveIterator;
    }
}

template <typename TInputImage, typename TOutputImage>
void
CombinedRelaxometryEstimationImageFilter <TInputImage,TOutputImage>
::SmoothB1Field()
{
    typedef itk::MultiplyImageFilter <OutputImageType,MaskImageType,OutputImageType> MultiplyFilterType;

    typename MultiplyFilterType::Pointer multiplier = MultiplyFilterType::New();
    multiplier->SetInput1(this->GetOutput(3));
    multiplier->SetInput2(this->GetComputationMask());
    multiplier->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

    multiplier->Update();

    typedef anima::SmoothingRecursiveYvvGaussianImageFilter <OutputImageType,OutputImageType> SmoothingFilterType;

    typename SmoothingFilterType::Pointer smoother = SmoothingFilterType::New();
    smoother->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    smoother->SetInput(multiplier->GetOutput());

    smoother->SetSigma(m_B1SmoothingSigma);

    smoother->Update();

    typedef anima::SmoothingRecursiveYvvGaussianImageFilter <MaskImageType,OutputImageType> SmoothingMaskFilterType;

    typename SmoothingMaskFilterType::Pointer maskSmoother = SmoothingMaskFilterType::New();
    maskSmoother->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    maskSmoother->SetInput(this->GetComputationMask());

    maskSmoother->SetSigma(m_B1SmoothingSigma);

    maskSmoother->Update();

    typename OutputImageType::Pointer tmpOutput = smoother->GetOutput();
    tmpOutput->DisconnectPipeline();

    typedef itk::ImageRegionIterator <OutputImageType> OutImageIteratorType;
    typedef itk::ImageRegionIterator <MaskImageType> MaskImageIteratorType;

    OutImageIteratorType smoothItr(tmpOutput,tmpOutput->GetLargestPossibleRegion());
    OutImageIteratorType maskSmoothItr(maskSmoother->GetOutput(),tmpOutput->GetLargestPossibleRegion());
    MaskImageIteratorType maskItr(this->GetComputationMask(),this->GetComputationMask()->GetLargestPossibleRegion());

    while (!maskItr.IsAtEnd())
    {
        if (maskItr.Get() == 0)
        {
            smoothItr.Set(1);

            ++smoothItr;
            ++maskItr;
            ++maskSmoothItr;
            continue;
        }

        double smoothVal = 1;
        if (maskSmoothItr.Get() != 0)
            smoothVal = smoothItr.Get() / maskSmoothItr.Get();

        smoothItr.Set(smoothVal);

        ++smoothItr;
        ++maskItr;
        ++maskSmoothItr;
    }

    this->SetNthOutput(3,tmpOutput);
}

} // end of namespace anima
