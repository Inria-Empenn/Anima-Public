#pragma once
#include "animaMultiT2RelaxometryEstimationImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

#include <animaNLOPTOptimizers.h>
#include <animaMultiT2EPGRelaxometryCostFunction.h>
#include <animaMeanAndVarianceImagesFilter.h>
#include <animaDekkerRootFindingAlgorithm.h>

namespace anima
{

template <class TPixelScalarType>
void
MultiT2RelaxometryEstimationImageFilter <TPixelScalarType>
::BeforeThreadedGenerateData()
{
    if ((m_RegularizationType == RegularizationType::NLTikhonov) && (!m_InitialT2Map || !m_InitialM0Map))
        itkExceptionMacro("Missing inputs for non-local estimation");

    Superclass::BeforeThreadedGenerateData();

    m_T2CompartmentValues.resize(m_NumberOfT2Compartments);
    double logStart = std::log(m_LowerT2Bound);
    double logEnd = std::log(m_UpperT2Bound);
    double step = (logEnd - logStart) / (m_NumberOfT2Compartments - 1.0);

    double logValue = logStart;
    m_T2CompartmentValues[0] = std::exp(logValue);
    for (unsigned int i = 1;i < m_NumberOfT2Compartments;++i)
    {
        logValue += step;
        m_T2CompartmentValues[i] = std::exp(logValue);
    }

    this->GetB1OutputImage()->FillBuffer(1.0);

    m_T2OutputImage = VectorOutputImageType::New();
    m_T2OutputImage->Initialize();
    m_T2OutputImage->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
    m_T2OutputImage->SetSpacing (this->GetInput(0)->GetSpacing());
    m_T2OutputImage->SetOrigin (this->GetInput(0)->GetOrigin());
    m_T2OutputImage->SetDirection (this->GetInput(0)->GetDirection());
    m_T2OutputImage->SetVectorLength(m_NumberOfT2Compartments);
    m_T2OutputImage->Allocate();

    OutputVectorType zero(m_NumberOfT2Compartments);
    zero.Fill(0.0);
    m_T2OutputImage->FillBuffer(zero);

    if (m_RegularizationType == RegularizationType::NLTikhonov)
    {
        this->PrepareNLPatchSearchers();

        typedef itk::ImageRegionConstIterator <InputImageType> InIteratorType;
        typedef itk::ImageRegionIterator <VectorOutputImageType> VectorIteratorType;

        InIteratorType  m0Iterator(m_InitialM0Map, m_InitialM0Map->GetLargestPossibleRegion());
        VectorIteratorType t2Iterator(m_InitialT2Map, m_InitialT2Map->GetLargestPossibleRegion());
        OutputVectorType tmpData;

        while (!m0Iterator.IsAtEnd())
        {
            tmpData = t2Iterator.Get();
            double m0Value = m0Iterator.Get();

            for (unsigned int i = 0;i < m_NumberOfT2Compartments;++i)
                tmpData[i] *= m0Value;

            t2Iterator.Set(tmpData);
            ++m0Iterator;
            ++t2Iterator;
        }
    }
}

template <class TPixelScalarType>
void
MultiT2RelaxometryEstimationImageFilter <TPixelScalarType>
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
MultiT2RelaxometryEstimationImageFilter <TPixelScalarType>
::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    typedef itk::ImageRegionConstIterator <InputImageType> ImageConstIteratorType;
    typedef itk::ImageRegionIterator <InputImageType> ImageIteratorType;
    typedef itk::ImageRegionIterator <VectorOutputImageType> VectorImageIteratorType;

    VectorImageIteratorType outT2Iterator(this->GetT2OutputImage(),outputRegionForThread);
    ImageIteratorType outM0Iterator(this->GetM0OutputImage(),outputRegionForThread);
    ImageIteratorType outMWFIterator(this->GetMWFOutputImage(),outputRegionForThread);
    ImageIteratorType outB1Iterator(this->GetB1OutputImage(),outputRegionForThread);
    ImageIteratorType outCostIterator(this->GetCostOutputImage(),outputRegionForThread);

    VectorImageIteratorType initT2Iterator;
    ImageIteratorType initM0Iterator;
    ImageIteratorType initB1Iterator;

    if (m_RegularizationType == RegularizationType::NLTikhonov)
    {
        initT2Iterator = VectorImageIteratorType(m_InitialT2Map, outputRegionForThread);
        initM0Iterator = ImageIteratorType(m_InitialM0Map, outputRegionForThread);
        initB1Iterator = ImageIteratorType(m_InitialB1Map, outputRegionForThread);
    }

    unsigned int numInputs = this->GetNumberOfIndexedInputs();
    std::vector <ImageConstIteratorType> inIterators;
    for (unsigned int i = 0;i < numInputs;++i)
        inIterators.push_back(ImageConstIteratorType(this->GetInput(i),outputRegionForThread));

    typedef itk::ImageRegionIterator <MaskImageType> MaskIteratorType;
    MaskIteratorType maskItr(this->GetComputationMask(),outputRegionForThread);

    ImageIteratorType t1MapItr;
    if (m_T1Map)
        t1MapItr = ImageIteratorType(m_T1Map,outputRegionForThread);

    itk::OptimizerParameters <double> signalValues(numInputs);
    itk::OptimizerParameters <double> signalValuesExtended(numInputs + m_NumberOfT2Compartments);
    itk::OptimizerParameters <double> t2OptimizedWeights(m_NumberOfT2Compartments);
    itk::OptimizerParameters <double> priorDistribution(m_NumberOfT2Compartments);
    priorDistribution.fill(0);

    OutputVectorType outputT2Weights(m_NumberOfT2Compartments);

    vnl_matrix <double> AMatrix(numInputs,m_NumberOfT2Compartments,0);
    vnl_matrix <double> AMatrixExtended(numInputs + m_NumberOfT2Compartments,m_NumberOfT2Compartments,0);

    typedef anima::NLOPTOptimizers B1OptimizerType;
    typedef anima::MultiT2EPGRelaxometryCostFunction B1CostFunctionType;

    typename B1CostFunctionType::Pointer cost = B1CostFunctionType::New();
    cost->SetEchoSpacing(m_EchoSpacing);
    cost->SetExcitationFlipAngle(m_T2ExcitationFlipAngle);

    unsigned int dimension = cost->GetNumberOfParameters();
    itk::Array<double> lowerBounds(dimension);
    itk::Array<double> upperBounds(dimension);
    B1OptimizerType::ParametersType p(dimension);
    lowerBounds[0] = 0.5 * m_T2FlipAngles[0];
    upperBounds[0] = 1.0 * m_T2FlipAngles[0];

    // NL specific variables
    std::vector <double> workDataWeights;
    std::vector <OutputVectorType> workDataSamples;

    unsigned int threadId = this->GetSafeThreadId();

    while (!maskItr.IsAtEnd())
    {
        outputT2Weights.Fill(0);

        if (maskItr.Get() == 0)
        {
            outT2Iterator.Set(outputT2Weights);
            outM0Iterator.Set(0);
            outMWFIterator.Set(0.0);
            outB1Iterator.Set(1.0);
            outCostIterator.Set(0.0);

            ++maskItr;
            ++outT2Iterator;
            ++outM0Iterator;
            ++outMWFIterator;
            ++outB1Iterator;
            ++outCostIterator;

            for (unsigned int i = 0;i < numInputs;++i)
                ++inIterators[i];

            if (m_T1Map)
                ++t1MapItr;

            if (m_RegularizationType == RegularizationType::NLTikhonov)
            {
                ++initB1Iterator;
                ++initT2Iterator;
                ++initM0Iterator;
            }

            continue;
        }

        double t1Value = 1000;
        double m0Value = 0.0;

        signalValuesExtended.fill(0.0);
        for (unsigned int i = 0;i < numInputs;++i)
        {
            signalValues[i] = inIterators[i].Get();
            signalValuesExtended[i] = signalValues[i];
        }

        if (m_T1Map)
        {
            t1Value = t1MapItr.Get();
            if (t1Value <= 0.0)
                t1Value = 1000;
        }

        cost->SetT1Value(t1Value);
        cost->SetT2Values(m_T2CompartmentValues);
        cost->SetT2RelaxometrySignals(signalValues);
        cost->SetUniformPulses(m_UniformPulses);
        if (!m_UniformPulses)
        {
            cost->SetPulseProfile(m_PulseProfile);
            cost->SetExcitationProfile(m_ExcitationProfile);
            cost->SetPixelWidth(m_PixelWidth);
        }

        double b1Value;

        if (m_InitialB1Map)
            b1Value = initB1Iterator.Get() * m_T2FlipAngles[0];
        else
            b1Value = 0.9 * m_T2FlipAngles[0];

        if (m_RegularizationType == RegularizationType::NLTikhonov)
        {
            outputT2Weights = initT2Iterator.Get();
            this->ComputeTikhonovPrior(maskItr.GetIndex(),outputT2Weights,m_NLPatchSearchers[threadId],priorDistribution,
                                       workDataWeights, workDataSamples);
        }

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

        b1Value = p[0] / m_T2FlipAngles[0];

        double residual = cost->GetValue(p);

        // Get data matrix back for regularization if needed
        AMatrix = cost->GetAMatrix();
        t2OptimizedWeights = cost->GetOptimizedT2Weights();
        m0Value = cost->GetOptimizedM0Value();

        AMatrixExtended.fill(0);
        for (unsigned int i = 0;i < m_NumberOfT2Compartments;++i)
        {
            for (unsigned int j = 0;j < numInputs;++j)
                AMatrixExtended(j,i) = AMatrix(j,i);
        }

        // Regularized NNLS with or without prior
        double normT2Weights = 0;
        for (unsigned int i = 0;i < m_NumberOfT2Compartments;++i)
            normT2Weights += (m0Value * t2OptimizedWeights[i] - priorDistribution[i]) * (m0Value * t2OptimizedWeights[i] - priorDistribution[i]);

        if (m_RegularizationType != RegularizationType::None)
        {
            RegularizationCostFunctionPointer regularizationCost = RegularizationCostFunctionType::New();
            regularizationCost->SetT2RelaxometrySignals(signalValuesExtended);
            regularizationCost->SetPriorDistribution(priorDistribution);
            regularizationCost->SetAMatrix(AMatrixExtended);
            regularizationCost->SetReferenceResidual(residual);
            regularizationCost->SetReferenceRatio(m_RegularizationRatio);
            regularizationCost->SetRegularizationType(m_RegularizationType);

            double lambda = std::sqrt(0.05 * residual / normT2Weights);
            residual = this->ComputeRegularizedSolution(regularizationCost,lambda,t2OptimizedWeights,m0Value);
        }

        outCostIterator.Set(residual);

        for (unsigned int i = 0;i < m_NumberOfT2Compartments;++i)
            outputT2Weights[i] = t2OptimizedWeights[i];

        outM0Iterator.Set(m0Value);
        outT2Iterator.Set(outputT2Weights);
        outCostIterator.Set(residual);
        double mwfValue = 0;
        for (unsigned int i = 0;i < m_NumberOfT2Compartments;++i)
        {
            if (m_T2CompartmentValues[i] > m_MyelinThreshold)
                break;

            mwfValue += outputT2Weights[i];
        }

        outMWFIterator.Set(mwfValue);
        outB1Iterator.Set(b1Value);

        this->IncrementNumberOfProcessedPoints();
        ++maskItr;
        ++outT2Iterator;
        ++outM0Iterator;
        ++outMWFIterator;
        ++outB1Iterator;
        ++outCostIterator;

        for (unsigned int i = 0;i < numInputs;++i)
            ++inIterators[i];

        if (m_T1Map)
            ++t1MapItr;

        if (m_RegularizationType == RegularizationType::NLTikhonov)
        {
            ++initB1Iterator;
            ++initT2Iterator;
            ++initM0Iterator;
        }
    }

    this->SafeReleaseThreadId(threadId);
}

template <class TPixelScalarType>
void
MultiT2RelaxometryEstimationImageFilter <TPixelScalarType>
::PrepareNLPatchSearchers()
{
    unsigned int numThreads = this->GetNumberOfWorkUnits();
    m_NLPatchSearchers.resize(numThreads);
    unsigned int maxAbsDisp = floor((double)(m_SearchNeighborhood / m_SearchStepSize)) * m_SearchStepSize;

    // Prepare mean and variance images
    unsigned int numInputs = this->GetNumberOfIndexedInputs();
    typedef anima::MeanAndVarianceImagesFilter<InputImageType, OutputImageType> MeanAndVarianceImagesFilterType;

    typename InputImageType::SizeType radius;
    for (unsigned int j = 0;j < InputImageType::ImageDimension;++j)
        radius[j] = m_PatchHalfSize;

    for (unsigned int i = 0;i < numInputs;++i)
    {
        typename MeanAndVarianceImagesFilterType::Pointer filter = MeanAndVarianceImagesFilterType::New();
        filter->SetInput(this->GetInput(i));

        filter->SetRadius(radius);
        filter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
        filter->Update();

        OutputImagePointer meanImage = filter->GetMeanImage();
        meanImage->DisconnectPipeline();
        OutputImagePointer varImage = filter->GetVarImage();
        varImage->DisconnectPipeline();

        for (unsigned int j = 0;j < numThreads;++j)
        {
            m_NLPatchSearchers[j].AddMeanImage(meanImage);
            m_NLPatchSearchers[j].AddVarImage(varImage);
            m_NLPatchSearchers[j].AddPatchTestImage(this->GetInput(i));
        }
    }

    // Prepare noise covariance values
    typename InputImageRegionType::IndexType baseIndex;
    typename InputImageRegionType::IndexType valueIndex;
    std::vector <double> noiseCovariances(numInputs,0);

    for (unsigned int i = 0;i < numInputs;++i)
    {
        typedef itk::ImageRegionConstIteratorWithIndex< InputImageType > InIteratorType;
        typename InputImageType::RegionType largestRegion = this->GetInput(i)->GetLargestPossibleRegion();
        InIteratorType  dataIterator (this->GetInput(i), largestRegion);

        double averageLocalSignal, diffSignal;
        double baseSignal;

        double averageCovariance = 0;
        unsigned int numEstimations = 0;
        unsigned int numLocalPixels = 2 * InputImageType::ImageDimension;

        while (!dataIterator.IsAtEnd())
        {
            baseSignal = static_cast<double>(dataIterator.Get());
            baseIndex = dataIterator.GetIndex();
            averageLocalSignal = 0;

            for (unsigned int d = 0; d < InputImageType::ImageDimension; ++d)
            {
                valueIndex = baseIndex;
                int tmpIndex = baseIndex[d] - m_LocalNeighborhood;
                valueIndex[d] = std::max(tmpIndex,0);
                averageLocalSignal += static_cast<double> (this->GetInput(i)->GetPixel(valueIndex));

                valueIndex = baseIndex;
                tmpIndex = baseIndex[d] + m_LocalNeighborhood;
                int maxIndex = largestRegion.GetSize()[d] - 1;
                valueIndex[d] = std::min(tmpIndex, maxIndex);
                averageLocalSignal += static_cast<double> (this->GetInput(i)->GetPixel(valueIndex));
            }

            averageLocalSignal /= numLocalPixels;
            diffSignal = sqrt(numLocalPixels / (numLocalPixels + 1.0)) * (baseSignal - averageLocalSignal);

            averageCovariance += diffSignal * diffSignal;

            ++numEstimations;
            ++dataIterator;
        }

        // Now divide by number of estimations and compute average variance
        noiseCovariances[i] = averageCovariance / numEstimations;
    }

    //Set searcher parameters
    for (unsigned int i = 0;i < numThreads;++i)
    {
        m_NLPatchSearchers[i].SetPatchHalfSize(m_PatchHalfSize);
        m_NLPatchSearchers[i].SetSearchStepSize(m_SearchStepSize);
        m_NLPatchSearchers[i].SetMaxAbsDisp(maxAbsDisp);
        m_NLPatchSearchers[i].SetInputImage(m_InitialT2Map);
        m_NLPatchSearchers[i].SetBetaParameter(m_BetaParameter);
        m_NLPatchSearchers[i].SetWeightThreshold(m_WeightThreshold);
        m_NLPatchSearchers[i].SetMeanMinThreshold(m_MeanMinThreshold);
        m_NLPatchSearchers[i].SetVarMinThreshold(m_VarMinThreshold);
        m_NLPatchSearchers[i].SetNoiseCovariances(noiseCovariances);
    }
}

template <class TPixelScalarType>
void
MultiT2RelaxometryEstimationImageFilter <TPixelScalarType>
::ComputeTikhonovPrior(const IndexType &refIndex, OutputVectorType &refDistribution,
                       PatchSearcherType &nlPatchSearcher, itk::OptimizerParameters <double> &priorDistribution,
                       std::vector <double> &workDataWeights, std::vector <OutputVectorType> &workDataSamples)
{
    nlPatchSearcher.UpdateAtPosition(refIndex);

    workDataWeights = nlPatchSearcher.GetDatabaseWeights();
    workDataSamples = nlPatchSearcher.GetDatabaseSamples();

    priorDistribution.fill(0);

    double sumWeights = 0;
    double maxWeight = -1;
    for (unsigned int i = 0;i < workDataSamples.size();++i)
    {
        bool isZero = true;
        for (unsigned int j = 0;j < m_NumberOfT2Compartments;++j)
        {
            if (workDataSamples[i][j] != 0)
            {
                isZero = false;
                break;
            }
        }

        if (isZero)
            continue;

        for (unsigned int j = 0;j < m_NumberOfT2Compartments;++j)
            priorDistribution[j] += workDataWeights[i] * workDataSamples[i][j];

        if (maxWeight < workDataWeights[i])
            maxWeight = workDataWeights[i];

        sumWeights += workDataWeights[i];
    }

    if (maxWeight < 0)
        maxWeight = 1.0;

    for (unsigned int j = 0;j < m_NumberOfT2Compartments;++j)
        priorDistribution[j] += maxWeight * refDistribution[j];

    sumWeights += maxWeight;
    priorDistribution /= sumWeights;
}

template <class TPixelScalarType>
double
MultiT2RelaxometryEstimationImageFilter <TPixelScalarType>
::ComputeRegularizedSolution(RegularizationCostFunctionType *regularizationCost, double &lambda,
                             itk::OptimizerParameters <double> &t2OptimizedWeights, double &m0Value)
{
    double lowerBound = 0.0;
    itk::OptimizerParameters <double> p(1);
    double zeroCost = 1.0 - m_RegularizationRatio;
    double upperBound = 0.0;

    // Find upper bound for lambda (and possibly lower bound)
    double ratio = -1.0;
    double oldRatio = ratio;
    bool ratioStall = false;
    while ((ratio < 0.0) && (!ratioStall))
    {
        oldRatio = ratio;
        p[0] = lambda;
        ratio = regularizationCost->GetValue(p);

        if (2.0 * std::abs(ratio - oldRatio) / std::abs(oldRatio + ratio) < 1.0e-2)
            ratioStall = true;

        if (ratio < 0.0)
        {
            lowerBound = lambda;
            lambda *= 2.0;
        }
    }

    if (ratio >= 0.0)
    {
        upperBound = lambda;

        DekkerRootFindingAlgorithm algorithm;

        algorithm.SetRootRelativeTolerance(1.0e-4);
        algorithm.SetCostFunctionTolerance(1.e-8);
        algorithm.SetRootFindingFunction(regularizationCost);
        algorithm.SetMaximumNumberOfIterations(100);
        algorithm.SetLowerBound(lowerBound);
        algorithm.SetUpperBound(upperBound);
        algorithm.SetFunctionValueAtInitialLowerBound(zeroCost);

        lambda = algorithm.Optimize();
        p[0] = lambda;
        regularizationCost->GetValue(p);

        t2OptimizedWeights = regularizationCost->GetOptimizedT2Weights();
        m0Value = regularizationCost->GetOptimizedM0Value();
    }

    return regularizationCost->GetCurrentResidual();
}

} // end namespace anima
