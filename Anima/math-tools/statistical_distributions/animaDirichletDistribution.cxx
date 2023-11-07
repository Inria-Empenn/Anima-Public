#include "animaDirichletDistribution.h"
#include <animaGammaFunctions.h>

#include <cmath>

#include <itkMacro.h>

namespace anima
{

bool DirichletDistribution::BelongsToSupport(const ValueType &x)
{
    unsigned int numParameters = x.size();

    double sumValue = 0.0;
    for (unsigned int i = 0;i < numParameters;++i)
    {
        double tmpValue = x[i];
        if (tmpValue < this->GetEpsilon() || tmpValue > 1.0 - this->GetEpsilon())
            return false;
        sumValue += tmpValue;
    }

    if (std::abs(sumValue - 1.0) > this->GetEpsilon())
        return false;
    
    return true;
}

void DirichletDistribution::SetConcentrationParameters(const std::vector<double> &val)
{
    unsigned int numParameters = val.size();
    
    for (unsigned int i = 0;i < numParameters;++i)
    {
        if (val[i] < this->GetEpsilon())
            throw itk::ExceptionObject(__FILE__, __LINE__, "The concentration parameters of a statistical distribution should be strictly positive.", ITK_LOCATION);
    }

    m_ConcentrationParameters = val;
}

double DirichletDistribution::GetDensity(const ValueType &x)
{
    if (!this->BelongsToSupport(x))
        return 0.0;

    return std::exp(this->GetLogDensity(x));
}

double DirichletDistribution::GetLogDensity(const ValueType &x)
{
    if (!this->BelongsToSupport(x))
        throw itk::ExceptionObject(__FILE__, __LINE__, "The log-density of the Dirichlet distribution is not defined for elements outside the simplex.", ITK_LOCATION);
    
    unsigned int numParameters = x.size();
    
    if (m_ConcentrationParameters.size() != numParameters)
        throw itk::ExceptionObject(__FILE__, __LINE__, "The input argument does not belong to the same simplex as the one on which the Dirichlet distribution is defined.", ITK_LOCATION);

    double logBeta = 0.0;
    double alphaZero = 0.0;
    double resValue = 0.0;
    for (unsigned int i = 0;i < numParameters;++i)
    {
        logBeta += std::lgamma(m_ConcentrationParameters[i]);
        alphaZero += m_ConcentrationParameters[i];
        resValue += (m_ConcentrationParameters[i] - 1.0) * std::log(x[i]);
    }
    logBeta -= std::lgamma(alphaZero);
    resValue -= logBeta;

    return resValue;
}

void DirichletDistribution::Fit(const SampleType &sample, const std::string &method)
{
    unsigned int numObservations = sample.size();
    unsigned int numParameters = sample[0].size();

    std::vector<double> alphaParameters(numParameters, 1.0);

    std::vector<bool> usefulValues(numObservations, false);
    ValueType sampleValue(numParameters);
    unsigned int numUsefulValues = 0;
    for (unsigned int i = 0;i < numObservations;++i)
    {
        for (unsigned int j = 0;j < numParameters;++j)
            sampleValue[j] = sample[i][j];
        if (this->BelongsToSupport(sampleValue))
        {
            usefulValues[i] = true;
            numUsefulValues++;
        }
    }

    if (numUsefulValues < 2)
    {
        this->SetConcentrationParameters(alphaParameters);
        return;
    }

    std::vector<double> logParameters(numParameters, 0.0);
    for (unsigned int j = 0;j < numParameters;++j)
    {
        for (unsigned int i = 0;i < numObservations;++i)
        {
            if (!usefulValues[i])
                continue;
            logParameters[j] += std::log(sample[i][j]);
        }

        logParameters[j] /= static_cast<double>(numUsefulValues);
    }

    // Compute initial sum of alpha's
    // From https://tminka.github.io/papers/dirichlet/minka-dirichlet.pdf
    // there are several options. We pick the one maximizing the first term in Eq. (4)
    
    // First compute empirical moments of order 1 and 2
    std::vector<double> firstMomentValues(numParameters);
    std::vector<double> secondMomentValues(numParameters);
    for (unsigned int j = 0;j < numParameters;++j)
    {
        double firstMoment = 0.0;
        double secondMoment = 0.0;
        for (unsigned int i = 0;i < numObservations;++i)
        {
            if (!usefulValues[i])
                continue;
            double tmpValue = sample[i][j];
            firstMoment += tmpValue;
            secondMoment += tmpValue * tmpValue;
        }
        firstMoment /= static_cast<double>(numUsefulValues);
        secondMoment /= static_cast<double>(numUsefulValues);

        firstMomentValues[j] = firstMoment;
        secondMomentValues[j] = secondMoment;
    }

    double alphaSum = 0.0;

    // Then, try Eq. (23) in https://tminka.github.io/papers/dirichlet/minka-dirichlet.pdf
    for (unsigned int j = 0;j < numParameters;++j)
    {
        double sumValue = 0.0;
        bool skipCandidate = false;
        for (unsigned int k = 0;k < numParameters;++k)
        {
            if (j == k)
                continue;
            
            double firstMoment = firstMomentValues[k];
            double secondMoment = secondMomentValues[k];
            
            double sampleVariance = secondMoment - firstMoment * firstMoment;
            if (sampleVariance < this->GetEpsilon())
            {
                skipCandidate = true;
                break;
            }

            double inLogValue = firstMoment * (1.0 - firstMoment) / sampleVariance - 1.0;
            if (inLogValue < this->GetEpsilon())
            {
                skipCandidate = true;
                break;
            }
            
            sumValue += std::log(inLogValue);
        }

        if (skipCandidate)
            continue;
        
        double candidateAlphaSum = std::exp(sumValue / (numParameters - 1.0));

        if (candidateAlphaSum > alphaSum)
            alphaSum = candidateAlphaSum;
    }

    // Then, try Eq. (21) from https://tminka.github.io/papers/dirichlet/minka-dirichlet.pdf
    for (unsigned int j = 0;j < numParameters;++j)
    {
        double firstMoment = firstMomentValues[j];
        double secondMoment = secondMomentValues[j];
        
        double sampleVariance = secondMoment - firstMoment * firstMoment;
        if (sampleVariance < this->GetEpsilon())
            continue;

        double candidateAlphaSum = (firstMoment - secondMoment) / sampleVariance;
        if (candidateAlphaSum < this->GetEpsilon())
            continue;
        
        if (candidateAlphaSum > alphaSum)
            alphaSum = candidateAlphaSum;
    }

    // Finally, if no solutions, output flat Dirichlet distribution
    if (alphaSum < this->GetEpsilon())
    {
        this->SetConcentrationParameters(alphaParameters);
        return;
    }
    
    double digammaValue = digamma(alphaSum);
    
    // Fixed point iteration method
    // Eq. (9) in https://tminka.github.io/papers/dirichlet/minka-dirichlet.pdf
    bool continueLoop = true;
    while (continueLoop)
    {
        double oldDigammaValue = digammaValue;

        for (unsigned int i = 0;i < numParameters;++i)
        {
            double psiNew = digammaValue + logParameters[i];
            alphaParameters[i] = anima::inverse_digamma(psiNew);
        }

        alphaSum = std::accumulate(alphaParameters.begin(), alphaParameters.end(), 0.0);
        digammaValue = digamma(alphaSum);
        continueLoop = std::abs(digammaValue - oldDigammaValue) > std::sqrt(std::numeric_limits<double>::epsilon());
    }

    m_TotalConcentration = alphaSum;
    m_MeanValues.resize(numParameters);
    for (unsigned int i = 0;i < numParameters;++i)
        m_MeanValues[i] = alphaParameters[i] / alphaSum;

    this->SetConcentrationParameters(alphaParameters);
}

void DirichletDistribution::Random(SampleType &sample, GeneratorType &generator)
{
    unsigned int numObservations = sample.size();
    unsigned int numParameters = sample[0].size();

    if (m_ConcentrationParameters.size() != numParameters)
        throw itk::ExceptionObject(__FILE__, __LINE__, "The requested sample does not belong to the same simplex as the one on which the Dirichlet distribution is defined.", ITK_LOCATION);

    BetaDistributionType betaDist;
    UniformDistributionType unifDist(0.0, 1.0);
    for (unsigned int i = 0;i < numObservations;++i)
    {
        double samplePartialSum = 0.0;
        for (unsigned int j = 0;j < numParameters - 1;++j)
        {
            double alphaPartialSum = 0.0;
            for (unsigned int k = j + 1;k < numParameters;++k)
                alphaPartialSum += m_ConcentrationParameters[k];

            betaDist = BetaDistributionType(m_ConcentrationParameters[j], alphaPartialSum);
            double phiValue = boost::math::quantile(betaDist, unifDist(generator));

            sample[i][j] = (1.0 - samplePartialSum) * phiValue;
            samplePartialSum += sample[i][j];
        }
        sample[i][numParameters - 1] = 1.0 - samplePartialSum;
    }
}

vnl_matrix<double> DirichletDistribution::GetCovarianceMatrix()
{
    unsigned int numComponents = m_ConcentrationParameters.size();
    vnl_matrix<double> covarianceMatrix(numComponents, numComponents);

    for (unsigned int i = 0;i < numComponents;++i)
    {
        covarianceMatrix(i, i) = m_MeanValues[i] * (1.0 - m_MeanValues[i]) / (1.0 + m_TotalConcentration);
        
        for (unsigned int j = i + 1;j < numComponents;++j)
        {
            double tmpValue = -1.0 * m_MeanValues[i] * m_MeanValues[j] / (1.0 + m_TotalConcentration);
            covarianceMatrix(i, j) = tmpValue;
            covarianceMatrix(j, i) = tmpValue;
        }
    }

    return covarianceMatrix;
}
    
} // end of namespace anima
