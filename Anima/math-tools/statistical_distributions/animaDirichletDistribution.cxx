#include "animaDirichletDistribution.h"
#include <cmath>
#include <limits>
#include <itkMacro.h>

namespace anima
{

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

    if (std::abs(sumValue - 1.0) < std::numeric_limits<double>::epsilon())
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

    if (std::abs(sumValue - 1.0) < std::numeric_limits<double>::epsilon())
        throw itk::ExceptionObject(__FILE__, __LINE__, "The log-density of the Dirichlet distribution is not defined for elements outside the simplex.", ITK_LOCATION);

    SingleValueType alphaParameters = this->GetConcentrationParameters();
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
    
}

void DirichletDistribution::Random(MultipleValueType &sample, GeneratorType &generator)
{
    unsigned int numObservations = sample.rows();
    unsigned int numParameters = sample.cols();

    SingleValueType alphaParameters = this->GetConcentrationParameters();
    if (alphaParameters.size() != numParameters)
        throw itk::ExceptionObject(__FILE__, __LINE__, "The requested sample does not belong to the same simplex as the one on which the Dirichlet distribution is defined.", ITK_LOCATION);

    BetaDistributionType betaDist;
    UniformDistributionType unifDist(0.0, 1.0);
    double samplePartialSum;
    for (unsigned int i = 0;i < numObservations;++i)
    {
        samplePartialSum = 0.0;
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
