#pragma once
#include "animaGammaMixtureT2RelaxometryEstimationImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

#include <animaNLOPTOptimizers.h>
#include <animaGammaMixtureT2RelaxometryCostFunction.h>

namespace anima
{

template <class TPixelScalarType>
void
GammaMixtureT2RelaxometryEstimationImageFilter <TPixelScalarType>
::BeforeThreadedGenerateData()
{
    Superclass::BeforeThreadedGenerateData();

    this->GetOutput(2)->FillBuffer(1.0);

    m_WeightsImage = VectorOutputImageType::New();
    m_WeightsImage->Initialize();
    m_WeightsImage->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
    m_WeightsImage->SetSpacing (this->GetInput(0)->GetSpacing());
    m_WeightsImage->SetOrigin (this->GetInput(0)->GetOrigin());
    m_WeightsImage->SetDirection (this->GetInput(0)->GetDirection());
    m_WeightsImage->SetVectorLength(3);
    m_WeightsImage->Allocate();

    OutputVectorType zero(3);
    zero.Fill(0.0);
    m_WeightsImage->FillBuffer(zero);

    m_MeanParamImage = VectorOutputImageType::New();
    m_MeanParamImage->Initialize();
    m_MeanParamImage->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
    m_MeanParamImage->SetSpacing (this->GetInput(0)->GetSpacing());
    m_MeanParamImage->SetOrigin (this->GetInput(0)->GetOrigin());
    m_MeanParamImage->SetDirection (this->GetInput(0)->GetDirection());
    m_MeanParamImage->SetVectorLength(3);
    m_MeanParamImage->Allocate();

    m_MeanParamImage->FillBuffer(zero);
}

template <class TPixelScalarType>
void
GammaMixtureT2RelaxometryEstimationImageFilter <TPixelScalarType>
::CheckComputationMask()
{
    if (this->GetComputationMask())
        return;

    typedef itk::ImageRegionConstIterator <InputImageType> IteratorType;
    typedef itk::ImageRegionIterator <MaskImageType> MaskIteratorType;

    std::vector <IteratorType> inItrs(this->GetNumberOfIndexedInputs());
    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
        inItrs[i] = IteratorType(this->GetInput(i),this->GetOutput(0)->GetLargestPossibleRegion());

    typename MaskImageType::Pointer maskImage = MaskImageType::New();
    maskImage->Initialize();
    maskImage->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
    maskImage->SetSpacing (this->GetInput(0)->GetSpacing());
    maskImage->SetOrigin (this->GetInput(0)->GetOrigin());
    maskImage->SetDirection (this->GetInput(0)->GetDirection());
    maskImage->Allocate();

    MaskIteratorType maskItr (maskImage,this->GetOutput(0)->GetLargestPossibleRegion());
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

template <class TPixelScalarType>
void
GammaMixtureT2RelaxometryEstimationImageFilter <TPixelScalarType>
::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    typedef itk::ImageRegionConstIterator <InputImageType> ImageConstIteratorType;
    typedef itk::ImageRegionIterator <InputImageType> ImageIteratorType;
    typedef itk::ImageRegionIterator <VectorOutputImageType> VectorImageIteratorType;

    VectorImageIteratorType outWeightsIterator(this->GetWeightsImage(),outputRegionForThread);
    VectorImageIteratorType outMeanParamsIterator(this->GetMeanParamImage(),outputRegionForThread);

    ImageIteratorType outM0Iterator(this->GetM0OutputImage(),outputRegionForThread);
    ImageIteratorType outMWFIterator(this->GetMWFOutputImage(),outputRegionForThread);
    ImageIteratorType outB1Iterator(this->GetB1OutputImage(),outputRegionForThread);
    ImageIteratorType outSigmaSqIterator(this->GetSigmaSquareOutputImage(),outputRegionForThread);

    unsigned int numInputs = this->GetNumberOfIndexedInputs();
    std::vector <ImageConstIteratorType> inIterators;
    for (unsigned int i = 0;i < numInputs;++i)
        inIterators.push_back(ImageConstIteratorType(this->GetInput(i),outputRegionForThread));

    typedef itk::ImageRegionIterator <MaskImageType> MaskIteratorType;
    MaskIteratorType maskItr(this->GetComputationMask(),outputRegionForThread);

    ImageIteratorType t1MapItr;
    if (m_T1Map)
        t1MapItr = ImageIteratorType(m_T1Map,outputRegionForThread);

    std::vector <double> VarParamsFixed(3, m_ShortT2Var);
    VarParamsFixed[1] = m_MediumT2Var;
    VarParamsFixed[2] = m_HighT2Var;

    std::vector <double> MeanParamsInit(3, m_ShortT2Mean);
    MeanParamsInit[1] = 110.0;
    MeanParamsInit[2] = m_HighT2Mean;

    using OptimizerType = anima::NLOPTOptimizers;
    using CostFunctionType = anima::B1GammaMixtureT2RelaxometryCostFunction;

    unsigned int n_compartments = 3;
    OptimizerType::ParametersType signalValues(numInputs);
    OptimizerType::ParametersType t2OptimizedWeights(n_compartments);

    OutputVectorType outputT2Weights(n_compartments);
    OutputVectorType outputMeanParams(n_compartments);

    typename CostFunctionType::Pointer cost = CostFunctionType::New();
    cost->SetEchoSpacing(m_EchoSpacing);
    cost->SetExcitationFlipAngle(m_T2ExcitationFlipAngle);
    cost->SetGammaMeans(MeanParamsInit);
    cost->SetGammaVariances(VarParamsFixed);
    cost->SetConstrainedParameters(m_ConstrainedParameters);
    cost->SetGammaIntegralTolerance(m_GammaIntegralTolerance);

    unsigned int dimension = cost->GetNumberOfParameters();

    itk::Array<double> lowerBounds(dimension);
    itk::Array<double> upperBounds(dimension);
    OptimizerType::ParametersType p(dimension);
    lowerBounds[0] = 1.0 * m_T2FlipAngles[0];
    upperBounds[0] = 2.0 * m_T2FlipAngles[0];

    // Initializations for L & U bounds and initial guesses.
    if (!m_ConstrainedParameters)
    {
        lowerBounds[1] = m_LowerShortT2;
        lowerBounds[2] = m_LowerMediumT2;
        lowerBounds[3] = m_LowerHighT2;

        upperBounds[1] = m_UpperShortT2;
        upperBounds[2] = m_UpperMediumT2;
        upperBounds[3] = m_UpperHighT2;
    }
    else
    {
        lowerBounds[1] = m_LowerMediumT2;
        upperBounds[1] = m_UpperMediumT2;
    }

    while (!maskItr.IsAtEnd())
    {
        outputT2Weights.Fill(0);
        outputMeanParams.Fill(0);

        if (maskItr.Get() == 0)
        {
            outWeightsIterator.Set(outputT2Weights);
            outMeanParamsIterator.Set(outputMeanParams);
            outM0Iterator.Set(0);
            outMWFIterator.Set(0.0);
            outB1Iterator.Set(0.0);
            outSigmaSqIterator.Set(0.0);

            ++maskItr;
            ++outWeightsIterator;
            ++outMeanParamsIterator;
            ++outM0Iterator;
            ++outMWFIterator;
            ++outB1Iterator;
            ++outSigmaSqIterator;

            for (unsigned int i = 0;i < numInputs;++i)
                ++inIterators[i];

            if (m_T1Map)
                ++t1MapItr;

            continue;
        }

        double t1Value = 1000.0;
        double m0Value = 0.0;

        for (unsigned int i = 0;i < numInputs;++i)
            signalValues[i] = inIterators[i].Get();

        if (m_T1Map)
        {
            t1Value = t1MapItr.Get();
            if (t1Value <= 0.0)
                t1Value = 1000.0;
        }

        cost->SetT1Value(t1Value);
        cost->SetT2RelaxometrySignals(signalValues);

        OptimizerType::Pointer opt = OptimizerType::New();
        opt->SetAlgorithm(NLOPT_LD_CCSAQ);
        opt->SetXTolRel(1.0e-5);
        opt->SetFTolRel(1.0e-7);
        opt->SetMaxEval(500);
        opt->SetVectorStorageSize(2000);

        opt->SetLowerBoundParameters(lowerBounds);
        opt->SetUpperBoundParameters(upperBounds);

        p[0] = 1.2 * m_T2FlipAngles[0];
        if (!m_ConstrainedParameters)
        {
            // Set Initial guesses
            p[1] = 20;
            p[2] = 110;
            p[3] = 2000;
        }
        else
            p[1] = 110;

        opt->SetInitialPosition(p);
        opt->SetMaximize(false);
        opt->SetCostFunction(cost);

        opt->StartOptimization();
        p = opt->GetCurrentPosition();

        t2OptimizedWeights = cost->GetOptimalT2Weights();

        // Get M0
        m0Value = 0;
        for (unsigned int i = 0;i < n_compartments;++i)
            m0Value += t2OptimizedWeights[i];

        for (unsigned int i = 0;i < n_compartments;++i)
        {
            if (m0Value > 0.0)
                outputT2Weights[i] = t2OptimizedWeights[i] / m0Value;
            else
                outputT2Weights[i] = 0.0;
        }

        outM0Iterator.Set(m0Value);
        outWeightsIterator.Set(outputT2Weights);

        if (!m_ConstrainedParameters)
        {
            for (unsigned int i = 0;i < 3;++i)
                outputMeanParams[i] = p[i + 1];
        }
        else
        {
            outputMeanParams[0] = m_ShortT2Mean;
            outputMeanParams[1] = p[1];
            outputMeanParams[2] = m_HighT2Mean;
        }

        outMeanParamsIterator.Set(outputMeanParams);

        double mwfValue = outputT2Weights[0];

        outMWFIterator.Set(mwfValue);

        double b1Value = p[0] / m_T2FlipAngles[0];
        outB1Iterator.Set(b1Value);

        outSigmaSqIterator.Set(cost->GetSigmaSquare());

        this->IncrementNumberOfProcessedPoints();
        ++maskItr;
        ++outWeightsIterator;
        ++outMeanParamsIterator;
        ++outM0Iterator;
        ++outMWFIterator;
        ++outB1Iterator;
        ++outSigmaSqIterator;

        for (unsigned int i = 0;i < numInputs;++i)
            ++inIterators[i];

        if (m_T1Map)
            ++t1MapItr;
    }
}

} // end namespace anima
