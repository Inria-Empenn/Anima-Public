#include "animaDirichletDistribution.h"
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

std::vector<std::vector<double>> CalculateMatrixLogarithm(const std::vector<std::vector<double>>& matrix) {
    std::vector<std::vector<double>> result;
    result.reserve(matrix.size());

    for (const auto& row : matrix) {
        std::vector<double> logRow;
        logRow.reserve(row.size());

        for (const auto& element : row) {
            if (element != 0)
                logRow.push_back(std::log(element));
            else
                logRow.push_back(0);
        }

        result.push_back(logRow);
    }

    return result;
}

void DirichletDistribution::Fit(const MultipleValueType &sample, const std::string &method)
{
    unsigned int numParameters = sample.cols();
    SingleValueType alphaParameters(numParameters, 1.0);
    unsigned int numObservations = sample.rows();
    unsigned int k;

    // calcul des log de la matrice
    std::vector<std::vector<double>> samplelog = calculateLogarithm(sample);
    
    std::vector<double> parametersLog;
    // Extract the columns into a separate range
    for (unsigned int i=0 ; i<numParameters, i++){
        size_t columnToSum = 1;
        std::vector<double> columnValues;
    
        for (const auto& row : matrix) {
            if (columnToSum < row.size()) {
            columnValues.push_back(row[columnToSum]);
            } else {
                // Handle column index out of range error
                throw std::out_of_range("Column index out of range");
            }
    }

        // Sum the columns values 
        double pk = std::accumulate(columnValues.begin(), columnValues.end(), 0.0);
        parametersLog.push_back(pk);
    } 

    double psi_new = 0;

    //boucle sur k
    for (unsigned int i = 0;i < 6;++i)
    {
        psi_new = digamma(std::accumulate(alphaParameters.begin(), alphaParameters.end(), 0)) + 1/numObservations * pk;
        alphaParameters[k] = inverse_digamma(psi_new);
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
