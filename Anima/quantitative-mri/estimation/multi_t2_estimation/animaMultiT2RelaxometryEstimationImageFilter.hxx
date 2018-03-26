#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

#include <animaEPGSignalSimulator.h>
#include <animaNLOPTOptimizers.h>
#include <animaMultiT2EPGRelaxometryCostFunction.h>
#include <animaMeanAndVarianceImagesFilter.h>

namespace anima
{

template <class TPixelScalarType>
void
MultiT2RelaxometryEstimationImageFilter <TPixelScalarType>
::BeforeThreadedGenerateData()
{
    if (m_NLEstimation && (!m_InitialT2Map || !m_InitialM0Map))
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

    if (m_NLEstimation)
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
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
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

    if (m_NLEstimation)
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

    T2VectorType signalValues(numInputs);
    T2VectorType signalValuesExtended(numInputs + m_NumberOfT2Compartments);
    T2VectorType t2OptimizedWeights(m_NumberOfT2Compartments);
    T2VectorType priorDistribution(m_NumberOfT2Compartments);
    priorDistribution.fill(0);

    OutputVectorType outputT2Weights(m_NumberOfT2Compartments);

    DataMatrixType AMatrix(numInputs,m_NumberOfT2Compartments,0);
    DataMatrixType AMatrixExtended(numInputs + m_NumberOfT2Compartments,m_NumberOfT2Compartments,0);

    anima::EPGSignalSimulator epgSimulator;
    epgSimulator.SetEchoSpacing(m_EchoSpacing);
    epgSimulator.SetNumberOfEchoes(numInputs);
    epgSimulator.SetB1OnExcitationAngle(m_B1OnExcitationAngle);
    epgSimulator.SetExcitationFlipAngle(m_T2ExcitationFlipAngle);
    epgSimulator.SetFlipAngle(m_T2FlipAngles[0]);

    anima::EPGSignalSimulator::RealVectorType epgSignalValues(numInputs);

    NNLSOptimizerPointer nnlsOpt = NNLSOptimizerType::New();

    typedef anima::NLOPTOptimizers B1OptimizerType;
    typedef anima::MultiT2EPGRelaxometryCostFunction B1CostFunctionType;

    typename B1CostFunctionType::Pointer cost = B1CostFunctionType::New();
    cost->SetEchoSpacing(m_EchoSpacing);
    cost->SetExcitationFlipAngle(m_T2ExcitationFlipAngle);
    cost->SetT2FlipAngles(m_T2FlipAngles);
    cost->SetB1OnExcitationAngle(m_B1OnExcitationAngle);
    cost->SetM0Value(1.0);

    unsigned int dimension = cost->GetNumberOfParameters();
    itk::Array<double> lowerBounds(dimension);
    itk::Array<double> upperBounds(dimension);
    B1OptimizerType::ParametersType p(dimension);
    lowerBounds[0] = 0;
    upperBounds[0] = 1;

    // NL specific variables
    std::vector <double> workDataWeights;
    std::vector <OutputVectorType> workDataSamples;

    while (!maskItr.IsAtEnd())
    {
        outputT2Weights.Fill(0);

        if (maskItr.Get() == 0)
        {
            outT2Iterator.Set(outputT2Weights);
            outM0Iterator.Set(0);
            outMWFIterator.Set(0.0);
            outB1Iterator.Set(0.0);
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

            if (m_NLEstimation)
            {
                ++initB1Iterator;
                ++initT2Iterator;
                ++initM0Iterator;
            }

            continue;
        }

        double previousB1Value = -1.0;
        double t1Value = 1000;
        double m0Value = 0.0;

        for (unsigned int i = 0;i < numInputs;++i)
        {
            signalValues[i] = inIterators[i].Get();
            signalValuesExtended[i] = signalValues[i];
        }

        if (m_T1Map)
            t1Value = t1MapItr.Get();
        cost->SetT1Value(t1Value);
        cost->SetT2Values(m_T2CompartmentValues);
        cost->SetT2RelaxometrySignals(signalValues);

        double b1Value;

        if (m_InitialB1Map)
            b1Value = initB1Iterator.Get();
        else
            b1Value = outB1Iterator.Get();

        if (m_NLEstimation)
        {
            outputT2Weights = initT2Iterator.Get();
            this->ComputeTikhonovPrior(maskItr.GetIndex(),outputT2Weights,m_NLPatchSearchers[threadId],priorDistribution,
                                       workDataWeights, workDataSamples);
        }

        unsigned int numGlobalIterations = 0;
        double residual = 0;

        while ((std::abs(b1Value - previousB1Value) > m_B1Tolerance)&&(numGlobalIterations < 100))
        {
            ++numGlobalIterations;
            previousB1Value = b1Value;
            // T2 weights estimation
            AMatrix.fill(0);
            AMatrixExtended.fill(0);
            for (unsigned int i = 0;i < m_NumberOfT2Compartments;++i)
            {
                epgSignalValues = epgSimulator.GetValue(t1Value,m_T2CompartmentValues[i],b1Value,1.0);
                for (unsigned int j = 0;j < numInputs;++j)
                {
                    AMatrix(j,i) = epgSignalValues[j];
                    AMatrixExtended(j,i) = epgSignalValues[j];
                }
            }

            // Regular NNLS optimization
            nnlsOpt->SetDataMatrix(AMatrix);
            nnlsOpt->SetPoints(signalValues);

            nnlsOpt->StartOptimization();
            t2OptimizedWeights = nnlsOpt->GetCurrentPosition();
            residual = nnlsOpt->GetCurrentResidual();

            // Regularized NNLS with or without prior
            double normT2Weights = 0;
            for (unsigned int i = 0;i < m_NumberOfT2Compartments;++i)
                normT2Weights += (t2OptimizedWeights[i] - priorDistribution[i]) * (t2OptimizedWeights[i] - priorDistribution[i]);

            if (m_RegularizationIntensity > 1.0)
            {
                double lambdaSq = (m_RegularizationIntensity - 1.0) * residual / normT2Weights;
                residual = this->ComputeTikhonovRegularizedSolution(nnlsOpt,AMatrixExtended,signalValuesExtended,
                                                                    lambdaSq,priorDistribution,t2OptimizedWeights);
            }

            outCostIterator.Set(residual);

            m0Value = 0.0;
            for (unsigned int i = 0;i < m_NumberOfT2Compartments;++i)
                m0Value += t2OptimizedWeights[i];

            for (unsigned int i = 0;i < m_NumberOfT2Compartments;++i)
                outputT2Weights[i] = t2OptimizedWeights[i] / m0Value;

            // B1 value estimation
            cost->SetT2Weights(t2OptimizedWeights);

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
            b1Value = p[0];
        }

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

        if (m_NLEstimation)
        {
            ++initB1Iterator;
            ++initT2Iterator;
            ++initM0Iterator;
        }
    }
}

template <class TPixelScalarType>
void
MultiT2RelaxometryEstimationImageFilter <TPixelScalarType>
::PrepareNLPatchSearchers()
{
    unsigned int numThreads = this->GetNumberOfThreads();
    m_NLPatchSearchers.resize(numThreads);
    unsigned int maxAbsDisp = floor((float)(m_SearchNeighborhood / m_SearchStepSize)) * m_SearchStepSize;

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
        filter->SetNumberOfThreads(this->GetNumberOfThreads());
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
                       PatchSearcherType &nlPatchSearcher, T2VectorType &priorDistribution,
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
::ComputeTikhonovRegularizedSolution(anima::NNLSOptimizer *nnlsOpt, DataMatrixType &AMatrix,
                                     T2VectorType &signalValues, double lambdaSq,
                                     T2VectorType &priorDistribution, T2VectorType &t2OptimizedWeights)
{
    unsigned int rowSize = AMatrix.rows() - AMatrix.cols();
    unsigned int colSize = AMatrix.cols();

    double lambda = std::sqrt(lambdaSq);

    for (unsigned int i = 0;i < colSize;++i)
    {
        signalValues[rowSize + i] = lambda * priorDistribution[i];
        for (unsigned int j = 0;j < colSize;++j)
        {
            if (i != j)
                AMatrix(rowSize + i,j) = 0;
            else
                AMatrix(rowSize + i,j) = lambda;
        }
    }

    nnlsOpt->SetDataMatrix(AMatrix);
    nnlsOpt->SetPoints(signalValues);

    nnlsOpt->StartOptimization();

    t2OptimizedWeights = nnlsOpt->GetCurrentPosition();

    return nnlsOpt->GetCurrentResidual();
}

} // end namespace anima
