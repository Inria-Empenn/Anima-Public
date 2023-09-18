#include "animaDirichletDistribution.h"
#include <animaGammaFunctions.h>
#include <cmath>
#include <limits>
#include <itkMacro.h>

namespace anima
{

void DirichletDistribution::SetConcentrationParameter(const SingleValueType val)
{
    unsigned int numParameters = val.size();

    for (unsigned int i = 0;i < numParameters;++i)
    {
        if (val[i] < std::numeric_limits<double>::epsilon())
            throw itk::ExceptionObject(__FILE__, __LINE__, "The concentration parameters of a statistical distribution should be strictly positive.", ITK_LOCATION);
    }

    this->BaseDistribution::SetConcentrationParameter(val);
}

double DirichletDistribution::GetDensity(const SingleValueType &x)
{
    unsigned int numParameters = x.size();

    double sumValue = 0.0;
    for (unsigned int i = 0;i < numParameters;++i)
    {
        double tmpValue = x[i];
        if (tmpValue < 0 || tmpValue > 1)
            return 0.0;
        sumValue += tmpValue;
    }

    if (std::abs(sumValue - 1.0) > std::numeric_limits<double>::epsilon())
        return 0.0;

    return std::exp(this->GetLogDensity(x));
}

double DirichletDistribution::GetLogDensity(const SingleValueType &x)
{
    unsigned int numParameters = x.size();

    double sumValue = 0.0;
    for (unsigned int i = 0;i < numParameters;++i)
    {
        double tmpValue = x[i];
        if (tmpValue < 0 || tmpValue > 1)
            throw itk::ExceptionObject(__FILE__, __LINE__, "The log-density of the Dirichlet distribution is not defined for elements outside the simplex.", ITK_LOCATION);
        sumValue += tmpValue;
    }

    if (std::abs(sumValue - 1.0) > std::numeric_limits<double>::epsilon())
        throw itk::ExceptionObject(__FILE__, __LINE__, "The log-density of the Dirichlet distribution is not defined for elements outside the simplex.", ITK_LOCATION);

    SingleValueType alphaParameters = this->GetConcentrationParameter();
    if (alphaParameters.size() != numParameters)
        throw itk::ExceptionObject(__FILE__, __LINE__, "The input argument does not belong to the same simplex as the one on which the Dirichlet distribution is defined.", ITK_LOCATION);

    double logBeta = 0.0;
    double alphaZero = 0.0;
    double resValue = 0.0;
    for (unsigned int i = 0;i < numParameters;++i)
    {
        logBeta += std::lgamma(alphaParameters[i]);
        alphaZero += alphaParameters[i];
        resValue += (alphaParameters[i] - 1.0) * std::log(x[i]);
    }
    logBeta -= std::lgamma(alphaZero);
    resValue -= logBeta;

    return resValue;
}

void DirichletDistribution::Fit(const MultipleValueType &sample, const std::string &method)
{
    unsigned int numObservations = sample.rows();
    unsigned int numParameters = sample.cols();

    SingleValueType alphaParameters(numParameters, 1.0);

    bool allZeros = true;
    SingleValueType logParameters(numParameters, 0.0);
    for (unsigned int j = 0;j < numParameters;++j)
    {
        unsigned int sumCount = 0;
        for (unsigned int i = 0;i < numObservations;++i)
        {
            double sampleValue = sample(i, j);
            if (sampleValue < std::numeric_limits<double>::epsilon())
                continue;
            logParameters[j] += std::log(sampleValue);
            sumCount++;
        }

        if (sumCount != 0)
        {
            logParameters[j] /= sumCount;
            allZeros = false;
        }
    }

    if (allZeros)
    {
        this->SetConcentrationParameter(alphaParameters);
        return;
    }

    // Compute initial sum of alpha's
    // Eq. (23) in https://tminka.github.io/papers/dirichlet/minka-dirichlet.pdf
    double sumValue = 0.0;
    unsigned int sumLength = 0;
    for (unsigned int j = 0;j < numParameters - 1;++j)
    {
        double firstMoment = 0.0;
        double secondMoment = 0.0;
        unsigned int realNumObservations = 0;

        for (unsigned int i = 0;i < numObservations;++i)
        {
            double tmpValue = sample(i, j);
            if (tmpValue < std::numeric_limits<double>::epsilon())
                continue;
            firstMoment += tmpValue;
            secondMoment += tmpValue * tmpValue;
            realNumObservations++;
        }

        if (realNumObservations < 2)
            continue;
        
        firstMoment /= realNumObservations;
        secondMoment /= realNumObservations;

        double sampleVariance = secondMoment - firstMoment * firstMoment;

        if (sampleVariance < std::numeric_limits<double>::epsilon())
            continue;

        double inLogValue = firstMoment * (1.0 - firstMoment) / sampleVariance - 1.0;

        if (inLogValue < std::numeric_limits<double>::epsilon())
            continue;
        
        sumValue += std::log(inLogValue);
        sumLength++;
    }

    if (sumLength == 0)
    {
        // Try Eq. (21) from https://tminka.github.io/papers/dirichlet/minka-dirichlet.pdf
        for (unsigned int j = 0;j < numParameters;++j)
        {
            double firstMoment = 0.0;
            double secondMoment = 0.0;
            unsigned int realNumObservations = 0;

            for (unsigned int i = 0;i < numObservations;++i)
            {
                double tmpValue = sample(i, j);
                if (tmpValue < std::numeric_limits<double>::epsilon())
                    continue;
                firstMoment += tmpValue;
                secondMoment += tmpValue * tmpValue;
                realNumObservations++;
            }

            if (realNumObservations < 2)
                continue;
            
            firstMoment /= realNumObservations;
            secondMoment /= realNumObservations;

            double sampleVariance = secondMoment - firstMoment * firstMoment;

            if (sampleVariance < std::numeric_limits<double>::epsilon())
                continue;
            
            double tmpAlphaSum = (firstMoment - secondMoment) / sampleVariance;

            if (tmpAlphaSum < sumValue || j == 0)
                sumValue = tmpAlphaSum;
        }

        if (sumValue == 0)
            sumValue = numParameters;

        sumValue = std::log(sumValue);
    }
    else
        sumValue /= sumLength;
    
    double alphaSum = std::exp(sumValue);
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

    this->SetConcentrationParameter(alphaParameters);
}

void DirichletDistribution::Random(MultipleValueType &sample, GeneratorType &generator)
{
    unsigned int numObservations = sample.rows();
    unsigned int numParameters = sample.cols();

    SingleValueType alphaParameters = this->GetConcentrationParameter();
    if (alphaParameters.size() != numParameters)
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
                alphaPartialSum += alphaParameters[k];

            betaDist = BetaDistributionType(alphaParameters[j], alphaPartialSum);
            double phiValue = boost::math::quantile(betaDist, unifDist(generator));

            sample(i, j) = (1.0 - samplePartialSum) * phiValue;
            samplePartialSum += sample(i, j);
        }
        sample(i, numParameters - 1) = 1.0 - samplePartialSum;
    }
}
    
} // end of namespace anima
