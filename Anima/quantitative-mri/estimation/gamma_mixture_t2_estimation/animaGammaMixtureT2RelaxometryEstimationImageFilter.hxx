#pragma once
#include "animaGammaMixtureT2RelaxometryEstimationImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

#include <animaNLOPTOptimizers.h>

#include <animaEPGSignalSimulator.h>
#include <animaGammaMixtureT2RelaxometryCostFunction.h>
#include <animaB1T2RelaxometryDistributionCostFunction.h>

#include <boost/math/distributions/gamma.hpp>

#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_svd.h>

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
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionConstIterator <InputImageType> ImageConstIteratorType;
    typedef itk::ImageRegionIterator <InputImageType> ImageIteratorType;
    typedef itk::ImageRegionIterator <VectorOutputImageType> VectorImageIteratorType;

    VectorImageIteratorType outWeightsIterator(this->GetWeightsImage(),outputRegionForThread);
    VectorImageIteratorType outMeanParamsIterator(this->GetMeanParamImage(),outputRegionForThread);

    ImageIteratorType outM0Iterator(this->GetM0OutputImage(),outputRegionForThread);
    ImageIteratorType outMWFIterator(this->GetMWFOutputImage(),outputRegionForThread);
    ImageIteratorType outB1Iterator(this->GetB1OutputImage(),outputRegionForThread);

    unsigned int numInputs = this->GetNumberOfIndexedInputs();
    std::vector <ImageConstIteratorType> inIterators;
    for (unsigned int i = 0;i < numInputs;++i)
        inIterators.push_back(ImageConstIteratorType(this->GetInput(i),outputRegionForThread));

    typedef itk::ImageRegionIterator <MaskImageType> MaskIteratorType;
    MaskIteratorType maskItr(this->GetComputationMask(),outputRegionForThread);

    ImageIteratorType t1MapItr;
    if (m_T1Map)
        t1MapItr = ImageIteratorType(m_T1Map,outputRegionForThread);

    T2VectorType signalValues(numInputs);
    unsigned int n_compartments = 3;
    T2VectorType t2OptimizedWeights(n_compartments);

    OutputVectorType outputT2Weights(n_compartments);
    OutputVectorType outputMeanParams(n_compartments);

    anima::EPGSignalSimulator epgSimulator;
    epgSimulator.SetEchoSpacing(m_EchoSpacing);
    epgSimulator.SetNumberOfEchoes(numInputs);
    epgSimulator.SetB1OnExcitationAngle(m_B1OnExcitationAngle);
    epgSimulator.SetExcitationFlipAngle(m_T2ExcitationFlipAngle);
    epgSimulator.SetFlipAngle(m_T2FlipAngles[0]);

    std::vector <anima::EPGSignalSimulator::RealVectorType> epgSignalValues;

    typedef anima::NLOPTOptimizers B1OptimizerType;
    typedef anima::B1T2RelaxometryDistributionCostFunction B1CostFunctionType;

    typename B1CostFunctionType::Pointer cost = B1CostFunctionType::New();
    cost->SetEchoSpacing(m_EchoSpacing);
    cost->SetExcitationFlipAngle(m_T2ExcitationFlipAngle);
    cost->SetT2FlipAngles(m_T2FlipAngles);
    cost->SetB1OnExcitationAngle(m_B1OnExcitationAngle);
    cost->SetM0Value(1.0);

    cost->SetLowerT2Bound(m_LowerT2Bound);
    cost->SetUpperT2Bound(m_UpperT2Bound);
    cost->SetT2IntegrationStep(m_T2IntegrationStep);

    unsigned int dimension = cost->GetNumberOfParameters();
    itk::Array<double> lowerBounds(dimension);
    itk::Array<double> upperBounds(dimension);
    B1OptimizerType::ParametersType p(dimension);
    lowerBounds[0] = 1;
    upperBounds[0] = 2;

    std::vector < std::vector <unsigned int> > sampledGammaT2Correspondences;
    std::vector < std::vector <double> > sampledGammaValues;
    std::vector <double> t2WorkingValues;

    // Initializations for L & U bounds and initial guesses.
    anima::NLOPTOptimizers::Pointer opt = anima::NLOPTOptimizers::New();
    
    anima::GammaMixtureT2RelaxometryCostFunction::Pointer cost_param = anima::GammaMixtureT2RelaxometryCostFunction::New();
    cost_param->SetConstrainedParameters(m_ConstrainedParameters);
    cost_param->SetNEchoes(numInputs);

    unsigned int numParametersGamma = cost_param->GetNumberOfParameters();
    NLOPTOptimizers::ParametersType optParameters(numParametersGamma);

    NLOPTOptimizers::ParametersType lowerBounds_GamParam(numParametersGamma);
    NLOPTOptimizers::ParametersType upperBounds_GamParam(numParametersGamma);

    if (!m_ConstrainedParameters)
    {
        lowerBounds_GamParam[0] = m_LowerShortT2;
        lowerBounds_GamParam[1] = m_LowerMediumT2;
        lowerBounds_GamParam[2] = m_LowerHighT2;

        upperBounds_GamParam[0] = m_UpperShortT2;
        upperBounds_GamParam[1] = m_UpperMediumT2;
        upperBounds_GamParam[2] = m_UpperHighT2;

        // Set Initial guesses
        optParameters[0] = 20;
        optParameters[1] = 110;
        optParameters[2] = 2000;
    }
    else
    {
        lowerBounds_GamParam[0] = m_LowerMediumT2;
        upperBounds_GamParam[0] = m_UpperMediumT2;

        // Set Initial guesses
        optParameters[0] = 110.0;
    }

    std::vector <double> VarParamsFixed(3, m_ShortT2Var);
    VarParamsFixed[1] = m_MediumT2Var;
    VarParamsFixed[2] = m_HighT2Var;

    std::vector <double> MeanParamsFixed(2, m_ShortT2Mean);
    MeanParamsFixed[1] = m_HighT2Mean;

    std::vector <double> SignalValuesVec(numInputs);
    std::vector <double> GamParamsVec(numParametersGamma);
    vnl_matrix <double> Lambda_Pinv;

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

            ++maskItr;
            ++outWeightsIterator;
            ++outMeanParamsIterator;
            ++outM0Iterator;
            ++outMWFIterator;
            ++outB1Iterator;

            for (unsigned int i = 0;i < numInputs;++i)
                ++inIterators[i];

            if (m_T1Map)
                ++t1MapItr;

            continue;
        }

        double previousB1Value = -1.0;
        double t1Value = 1000;
        double m0Value = 1.0;

        for (unsigned int i = 0;i < numInputs;++i)
            signalValues[i] = inIterators[i].Get();

        if (m_T1Map)
            t1Value = t1MapItr.Get();

        cost->SetT1Value(t1Value);
        cost->SetT2RelaxometrySignals(signalValues);

        double b1Value = 1.2;
        unsigned int numGlobalIterations = 0;

        for (unsigned int i = 0; i < numInputs; ++i)
            SignalValuesVec[i] = signalValues[i];

        // Generate the t2-indices for Parameter Optimization
        unsigned int numSteps = std::ceil((m_UpperT2Bound - m_LowerT2Bound) / m_T2IntegrationStep) + 1;
        for (unsigned int i = 0; i < numSteps; ++i)
            t2WorkingValues.push_back(m_LowerT2Bound + i * m_T2IntegrationStep);

        double previousCostValue = 0;
        double costValue = 1;

        while ((!this->endConditionReached(b1Value,previousB1Value,costValue,previousCostValue))&&(numGlobalIterations < 100))
        {
            ++numGlobalIterations;
            previousB1Value = b1Value;
            previousCostValue = costValue;

            epgSignalValues.resize(numSteps);
            for (unsigned int i = 0; i < numSteps; ++i)
                epgSignalValues[i] = epgSimulator.GetValue(t1Value,t2WorkingValues[i],b1Value,1.0);

            cost_param->SetNumSteps(numSteps);
            cost_param->SetT2WorkingValues(t2WorkingValues);
            cost_param->SetEPGSignalValues(epgSignalValues);
            cost_param->SetSignalValues(SignalValuesVec);
            cost_param->SetMeanParam(MeanParamsFixed);
            cost_param->SetVarParam(VarParamsFixed);
            opt->SetCostFunction(cost_param);

            //Set the 'opt' parameters
            opt->SetAlgorithm(NLOPT_LD_CCSAQ);
            opt->SetXTolRel(1.0e-4);
            opt->SetFTolRel(1.0e-6);
            opt->SetMaxEval(500);
            opt->SetVectorStorageSize(2000);

            opt->SetLowerBoundParameters(lowerBounds_GamParam);
            opt->SetUpperBoundParameters(upperBounds_GamParam);
            opt->SetInitialPosition(optParameters);
            opt->SetMaximize(false);

            opt->StartOptimization();
            optParameters = opt->GetCurrentPosition();

            for (unsigned int i = 0;i < numParametersGamma;++i)
                GamParamsVec[i] = optParameters[i];

            costValue = cost_param->GetValue(optParameters);
            Lambda_Pinv = cost_param->GetLambda_pinv();

            t2OptimizedWeights.fill(0);

            vnl_vector<double> SignalValues_vnl_format(SignalValuesVec.size());
            for (unsigned int i = 0; i < SignalValuesVec.size(); ++i)
                SignalValues_vnl_format[i] = SignalValuesVec[i];

            t2OptimizedWeights = Lambda_Pinv * SignalValues_vnl_format;
            m0Value = 0;
            for (unsigned int i = 0; i < t2OptimizedWeights.size(); ++i)
                m0Value += t2OptimizedWeights[i];

            for (unsigned int i = 0; i < t2OptimizedWeights.size(); ++i)
            {
                t2OptimizedWeights[i] /= m0Value;
                outputT2Weights[i] = t2OptimizedWeights[i];
            }

            this->PrepareGammaValues(sampledGammaValues,sampledGammaT2Correspondences,t2WorkingValues,GamParamsVec,VarParamsFixed);
            cost->SetT2DistributionSamples(sampledGammaValues);
            cost->SetT2WorkingValues(t2WorkingValues);
            cost->SetDistributionSamplesT2Correspondences(sampledGammaT2Correspondences);

            cost->SetT2Weights(t2OptimizedWeights);
            cost->SetM0Value(m0Value);

            B1OptimizerType::Pointer b1Optimizer = B1OptimizerType::New();
            b1Optimizer->SetAlgorithm(NLOPT_LN_BOBYQA);
            b1Optimizer->SetXTolRel(1.0e-4);
            b1Optimizer->SetFTolRel(1.0e-6);
            b1Optimizer->SetMaxEval(500);
            b1Optimizer->SetVectorStorageSize(2000);

            b1Optimizer->SetLowerBoundParameters(lowerBounds);
            b1Optimizer->SetUpperBoundParameters(upperBounds);
            p[0] = b1Value;
            b1Optimizer->SetInitialPosition(p);
            b1Optimizer->SetMaximize(false);
            b1Optimizer->SetCostFunction(cost);

            b1Optimizer->StartOptimization();
            p = b1Optimizer->GetCurrentPosition();

            costValue = cost->GetValue(p);
            b1Value = p[0];
        }

        if (!m_ConstrainedParameters)
        {
            for (unsigned int i = 0; i < optParameters.size(); ++i)
                outputMeanParams[i] = optParameters[i];
        }
        else
        {
            outputMeanParams[0] = m_ShortT2Mean;
            outputMeanParams[1] = optParameters[0];
            outputMeanParams[2] = m_HighT2Mean;
        }

        outM0Iterator.Set(m0Value);
        outWeightsIterator.Set(outputT2Weights);
        outMeanParamsIterator.Set(outputMeanParams);

        double mwfValue = outputT2Weights[0];

        outMWFIterator.Set(mwfValue);
        outB1Iterator.Set(b1Value);

        ++maskItr;
        ++outWeightsIterator;
        ++outMeanParamsIterator;
        ++outM0Iterator;
        ++outMWFIterator;
        ++outB1Iterator;

        for (unsigned int i = 0;i < numInputs;++i)
            ++inIterators[i];

        if (m_T1Map)
            ++t1MapItr;
    }
}

template <class TPixelScalarType>
bool
GammaMixtureT2RelaxometryEstimationImageFilter <TPixelScalarType>
::endConditionReached(double b1Value, double previousB1Value, double costValue, double previousCostValue)
{
    if (std::abs(costValue - previousCostValue) / costValue < m_CostTolerance)
        return true;

    if (std::abs(b1Value - previousB1Value) > m_B1Tolerance)
        return false;

    return true;
}

// Gamma-Distribution Preparation: Prefixed Alpha and Lambda Values
template <class TPixelScalarType>
void
GammaMixtureT2RelaxometryEstimationImageFilter <TPixelScalarType>
::PrepareGammaValues(std::vector < std::vector <double> > &sampledGammaValues, std::vector < std::vector <unsigned int> > &sampledGammaT2Correspondences,
                     std::vector <double> &t2WorkingValues, std::vector <double> &optimizedGamParams, std::vector <double> &gammaVariance)
{
    // Allocation of the Gamma Distribution Parameters
    unsigned int numberOfGamma = 3;

    std::vector <double> gammaMean(numberOfGamma);

    if (m_ConstrainedParameters)
    {
        gammaMean[0] = m_ShortT2Mean;
        gammaMean[1] = optimizedGamParams[0];
        gammaMean[2] = m_HighT2Mean;
    }
    else
    {
        for (unsigned int i = 0;i < numberOfGamma;++i)
            gammaMean[i] = optimizedGamParams[i];
    }

    sampledGammaValues.resize(numberOfGamma);
    sampledGammaT2Correspondences.resize(numberOfGamma);
    t2WorkingValues.clear();

    unsigned int numSteps = std::ceil((m_UpperT2Bound - m_LowerT2Bound) / m_T2IntegrationStep) + 1;
    std::vector <unsigned int> tmpVec;
    std::vector < std::vector <unsigned int> > usedIndexes(numSteps,tmpVec);

    for (unsigned int i = 0;i < numSteps;++i)
        t2WorkingValues.push_back(m_LowerT2Bound + i * m_T2IntegrationStep);

    std::vector <double> gammaVector(numSteps,0.0);
    std::vector <unsigned int> gammaT2Indexes(numSteps,0);

    for (unsigned int j = 0;j < numberOfGamma;++j)
    {
        double shape = gammaMean[j] * gammaMean[j] / gammaVariance[j];
        double scale = gammaVariance[j] / gammaMean[j];

        for (unsigned int i = 0;i < numSteps;++i)
        {
            usedIndexes[i].push_back(j);
            double position = t2WorkingValues[i];
            double internalValue = boost::math::gamma_p_derivative(shape, position / scale) / scale;
            gammaVector[i] = internalValue;
            gammaT2Indexes[i] = i;
        }

        sampledGammaValues[j] = gammaVector;
        sampledGammaT2Correspondences[j] = gammaT2Indexes;
    }
}

} // end namespace anima
