#pragma once
#include "animaGMMT2RelaxometryEstimationImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

#include <animaEPGSignalSimulator.h>
#include <animaNLOPTOptimizers.h>
#include <animaB1T2RelaxometryDistributionCostFunction.h>

#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/erf.hpp>

namespace anima
{

template <class TPixelScalarType>
void
GMMT2RelaxometryEstimationImageFilter <TPixelScalarType>
::SetGaussianMeans(std::string fileName)
{
    std::ifstream fileIn(fileName.c_str());

    if (!fileIn.is_open())
    {
        std::string errorMsg = "Unable to read file: ";
        errorMsg += fileName;
        throw itk::ExceptionObject(__FILE__,__LINE__,errorMsg,ITK_LOCATION);
    }

    m_GaussianMeans.clear();

    while (!fileIn.eof())
    {
        char tmpStr[2048];
        fileIn.getline(tmpStr,2048);

        std::string workStr(tmpStr);
        if (workStr == "")
            continue;

        boost::algorithm::trim_right(workStr);
        std::istringstream iss(workStr);
        std::string shortStr;
        iss >> shortStr;
        m_GaussianMeans.push_back(boost::lexical_cast<double> (shortStr));
    }

    fileIn.close();
}

template <class TPixelScalarType>
void
GMMT2RelaxometryEstimationImageFilter <TPixelScalarType>
::SetGaussianVariances(std::string fileName)
{
    std::ifstream fileIn(fileName.c_str());

    if (!fileIn.is_open())
    {
        std::string errorMsg = "Unable to read file: ";
        errorMsg += fileName;
        throw itk::ExceptionObject(__FILE__,__LINE__,errorMsg,ITK_LOCATION);
    }

    m_GaussianVariances.clear();

    while (!fileIn.eof())
    {
        char tmpStr[2048];
        fileIn.getline(tmpStr,2048);

        std::string workStr(tmpStr);
        if (workStr == "")
            continue;

        boost::algorithm::trim_right(workStr);
        std::istringstream iss(workStr);
        std::string shortStr;
        iss >> shortStr;
        m_GaussianVariances.push_back(boost::lexical_cast<double> (shortStr));
    }

    fileIn.close();
}

template <class TPixelScalarType>
void
GMMT2RelaxometryEstimationImageFilter <TPixelScalarType>
::BeforeThreadedGenerateData()
{
    if (m_GaussianVariances.size() != 3)
        itkExceptionMacro("Number of variance parameters has to be 3...");

    Superclass::BeforeThreadedGenerateData();

    this->GetB1OutputImage()->FillBuffer(1.0);

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

    this->PrepareGaussianValues();
}

template <class TPixelScalarType>
void
GMMT2RelaxometryEstimationImageFilter <TPixelScalarType>
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
GMMT2RelaxometryEstimationImageFilter <TPixelScalarType>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionConstIterator <InputImageType> ImageConstIteratorType;
    typedef itk::ImageRegionIterator <InputImageType> ImageIteratorType;
    typedef itk::ImageRegionIterator <VectorOutputImageType> VectorImageIteratorType;

    VectorImageIteratorType outWeightsIterator(this->GetWeightsImage(),outputRegionForThread);
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
    unsigned int numberOfGaussians = m_GaussianMeans.size();
    T2VectorType t2OptimizedWeights(numberOfGaussians);

    OutputVectorType outputT2Weights(numberOfGaussians);
    DataMatrixType AMatrix(numInputs,numberOfGaussians,0);

    anima::EPGSignalSimulator epgSimulator;
    epgSimulator.SetEchoSpacing(m_EchoSpacing);
    epgSimulator.SetNumberOfEchoes(numInputs);
    epgSimulator.SetB1OnExcitationAngle(m_B1OnExcitationAngle);
    epgSimulator.SetExcitationFlipAngle(m_T2ExcitationFlipAngle);
    epgSimulator.SetFlipAngle(m_T2FlipAngles[0]);

    std::vector <anima::EPGSignalSimulator::RealVectorType> epgSignalValues;

    NNLSOptimizerPointer nnlsOpt = NNLSOptimizerType::New();

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

    cost->SetT2DistributionSamples(m_SampledGaussianValues);
    cost->SetT2WorkingValues(m_T2WorkingValues);
    cost->SetDistributionSamplesT2Correspondences(m_SampledGaussT2Correspondences);

    unsigned int numSteps = m_T2WorkingValues.size();
    epgSignalValues.resize(numSteps);

    while (!maskItr.IsAtEnd())
    {
        outputT2Weights.Fill(0);

        if (maskItr.Get() == 0)
        {
            outWeightsIterator.Set(outputT2Weights);
            outM0Iterator.Set(0);
            outMWFIterator.Set(0.0);
            outB1Iterator.Set(0.0);

            ++maskItr;
            ++outWeightsIterator;
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
        double m0Value = 0.0;

        for (unsigned int i = 0;i < numInputs;++i)
            signalValues[i] = inIterators[i].Get();

        if (m_T1Map)
            t1Value = t1MapItr.Get();
        cost->SetT1Value(t1Value);
        cost->SetT2RelaxometrySignals(signalValues);

        double b1Value = outB1Iterator.Get();
        unsigned int numGlobalIterations = 0;

        while ((!this->endConditionReached(b1Value,previousB1Value))&&(numGlobalIterations < 100))
        {
            ++numGlobalIterations;
            previousB1Value = b1Value;
            // T2 weights estimation
            AMatrix.fill(0);

            for (unsigned int i = 0;i < numSteps;++i)
                epgSignalValues[i] = epgSimulator.GetValue(t1Value,m_T2WorkingValues[i],b1Value,1.0);

            for (unsigned int i = 0;i < numberOfGaussians;++i)
            {
                unsigned int numSamples = m_SampledGaussianValues[i].size();
                for (unsigned int j = 0;j < numInputs;++j)
                {
                    double integralValue = m_SampledGaussianValues[i][0] * epgSignalValues[m_SampledGaussT2Correspondences[i][0]][j] * m_T2IntegrationStep / 2.0;

                    for (unsigned int k = 1;k < numSamples;++k)
                        integralValue += (m_SampledGaussianValues[i][k] * epgSignalValues[m_SampledGaussT2Correspondences[i][k]][j] + m_SampledGaussianValues[i][k-1] * epgSignalValues[m_SampledGaussT2Correspondences[i][k-1]][j]) * m_T2IntegrationStep / 2.0;

                    AMatrix(j,i) = integralValue;
                }
            }

            // Regular NNLS optimization
            nnlsOpt->SetDataMatrix(AMatrix);
            nnlsOpt->SetPoints(signalValues);

            nnlsOpt->StartOptimization();
            t2OptimizedWeights = nnlsOpt->GetCurrentPosition();

            m0Value = 0.0;
            for (unsigned int i = 0;i < numberOfGaussians;++i)
                m0Value += t2OptimizedWeights[i];

            for (unsigned int i = 0;i < numberOfGaussians;++i)
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
        outWeightsIterator.Set(outputT2Weights);

        double mwfValue = outputT2Weights[0];

        outMWFIterator.Set(mwfValue);
        outB1Iterator.Set(b1Value);

        ++maskItr;
        ++outWeightsIterator;
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
GMMT2RelaxometryEstimationImageFilter <TPixelScalarType>
::endConditionReached(double b1Value, double previousB1Value)
{
    if (std::abs(b1Value - previousB1Value) > m_B1Tolerance)
        return false;

    return true;
}

template <class TPixelScalarType>
void
GMMT2RelaxometryEstimationImageFilter <TPixelScalarType>
::PrepareGaussianValues()
{
    unsigned int numberOfGaussians = m_GaussianMeans.size();

    m_SampledGaussianValues.resize(numberOfGaussians);
    m_SampledGaussT2Correspondences.resize(numberOfGaussians);
    m_T2WorkingValues.clear();

    unsigned int maxNumSteps = std::ceil((m_UpperT2Bound - m_LowerT2Bound) / m_T2IntegrationStep) + 1;
    std::vector <unsigned int> tmpVec;
    std::vector < std::vector <unsigned int> > usedIndexes(maxNumSteps,tmpVec);

    double widthGaussianApproximation = boost::math::erf_inv(1.0 - m_GaussianIntegralTolerance) * std::sqrt(2.0);

    for (unsigned int j = 0;j < numberOfGaussians;++j)
    {
        double minValueGaussian = std::max(m_LowerT2Bound, m_GaussianMeans[j] - widthGaussianApproximation * std::sqrt(m_GaussianVariances[j]));
        double maxValueGaussian = std::min(m_UpperT2Bound, m_GaussianMeans[j] + widthGaussianApproximation * std::sqrt(m_GaussianVariances[j]));

        unsigned int iMin = std::ceil((minValueGaussian - m_LowerT2Bound) / m_T2IntegrationStep);
        unsigned int iMax = std::floor((maxValueGaussian - m_LowerT2Bound) / m_T2IntegrationStep);

        unsigned int numSteps = iMax - iMin + 1;
        std::vector <double> gaussianVector(numSteps,0.0);
        std::vector <unsigned int> gaussT2Indexes(numSteps,0);
        double vectorSum = 0;

        for (unsigned int i = iMin;i <= iMax;++i)
        {
            usedIndexes[i].push_back(j);
            double position = i * m_T2IntegrationStep + m_LowerT2Bound;
            double internalValue = (position - m_GaussianMeans[j]) * (position - m_GaussianMeans[j]) / (2.0 * m_GaussianVariances[j]);
            gaussianVector[i - iMin] = std::exp(- internalValue);
            gaussT2Indexes[i - iMin] = i;

            if (i > iMin)
                vectorSum += (gaussianVector[i - iMin] + gaussianVector[i - iMin - 1]) * m_T2IntegrationStep / 2.0;
            else
                vectorSum += gaussianVector[i - iMin] * m_T2IntegrationStep / 2.0;
        }

        for (unsigned int i = 0;i < numSteps;++i)
            gaussianVector[i] /= vectorSum;

        m_SampledGaussianValues[j] = gaussianVector;
        m_SampledGaussT2Correspondences[j] = gaussT2Indexes;
    }

    for (unsigned int i = 0;i < maxNumSteps;++i)
    {
        if (usedIndexes[i].size() != 0)
        {
            for (unsigned int j = 0;j < usedIndexes[i].size();++j)
            {
                double minValueGaussian = std::max(m_LowerT2Bound, m_GaussianMeans[usedIndexes[i][j]] - widthGaussianApproximation * std::sqrt(m_GaussianVariances[usedIndexes[i][j]]));
                unsigned int iMin = std::ceil((minValueGaussian - m_LowerT2Bound) / m_T2IntegrationStep);

                m_SampledGaussT2Correspondences[usedIndexes[i][j]][i - iMin] = m_T2WorkingValues.size();
            }

            m_T2WorkingValues.push_back(m_LowerT2Bound + i * m_T2IntegrationStep);
        }
    }
}

} // end namespace anima
