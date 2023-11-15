#include "animaMultivariateNormalDistribution.h"
#include <animaBaseTensorTools.h>

#include <itkMacro.h>

#include <vnl/algo/vnl_determinant.h>

namespace anima
{
    void MultivariateNormalDistribution::SetCovarianceMatrixParameter(const MatrixType &val)
    {
        unsigned int dimValue = val.rows();

        if (val.cols() != dimValue)
            throw itk::ExceptionObject(__FILE__, __LINE__, "The covariance matrix is not square.", ITK_LOCATION);

        bool isSymmetric = true;
        for (unsigned int i = 0; i < dimValue; ++i)
        {
            for (unsigned int j = i + 1; j < dimValue; ++j)
            {
                if (val.get(i, j) != val.get(j, i))
                    throw itk::ExceptionObject(__FILE__, __LINE__, "The covariance matrix is not symmetric.", ITK_LOCATION);
            }
        }

        m_CovarianceMatrixDeterminant = vnl_determinant<double>(val);
        if (m_CovarianceMatrixDeterminant < this->GetEpsilon())
            throw itk::ExceptionObject(__FILE__, __LINE__, "The covariance matrix is not definite positive.", ITK_LOCATION);

        m_CovarianceMatrixParameter = val;

        MatrixType eigenVecs(dimValue, dimValue);
        vnl_diag_matrix<double> eigenVals(dimValue), sqrtEigenVals(dimValue), invEigenVals(dimValue);

        itk::SymmetricEigenAnalysis<MatrixType, vnl_diag_matrix<double>, MatrixType> eigenComputer(dimValue);
        eigenComputer.ComputeEigenValuesAndVectors(m_CovarianceMatrixParameter, eigenVals, eigenVecs);

        for (unsigned int i = 0; i < dimValue; ++i)
        {
            sqrtEigenVals[i] = std::sqrt(eigenVals[i]);
            invEigenVals[i] = 1.0 / eigenVals[i];
        }

        anima::RecomposeTensor(sqrtEigenVals, eigenVecs, m_SqrtCovarianceMatrixParameter);
        anima::RecomposeTensor(invEigenVals, eigenVecs, m_PrecisionMatrixParameter);
    }

    double MultivariateNormalDistribution::GetDensity(const ValueType &x)
    {
        if (!this->BelongsToSupport(x))
            return 0.0;

        return std::exp(this->GetLogDensity(x));
    }

    double MultivariateNormalDistribution::GetLogDensity(const ValueType &x)
    {
        if (!this->BelongsToSupport(x))
            throw itk::ExceptionObject(__FILE__, __LINE__, "The log-density is not defined outside the support.", ITK_LOCATION);

        unsigned int numParameters = x.size();

        if (m_MeanParameter.size() != numParameters || m_CovarianceMatrixParameter.rows() != numParameters)
            throw itk::ExceptionObject(__FILE__, __LINE__, "The input argument does not belong to the same vector space as the one on which the normal distribution is defined.", ITK_LOCATION);

        double resValue = -static_cast<double>(numParameters) / 2.0 * std::log(2.0 * M_PI);
        resValue -= 0.5 * std::log(m_CovarianceMatrixDeterminant);

        double tmpValue = 0.0;
        for (unsigned int i = 0; i < numParameters; ++i)
            for (unsigned int j = 0; j < numParameters; ++j)
                tmpValue += m_PrecisionMatrixParameter.get(i, j) * (x[i] - m_MeanParameter[i]) * (x[j] * m_MeanParameter[j]);

        resValue -= 0.5 * tmpValue;

        return resValue;
    }

    double MultivariateNormalDistribution::GetCumulative(const ValueType &x)
    {
        /**
         * \fn double MultivariateNormalDistribution::GetCumulative(const std::vector<double> &x)
         *
         * \author Aymeric Stamm
         * \date November 2023
         *
         * \param x A numeric vector specifying a point in R^p.
         *
         * \return A numeric value storing the value of the CDF of the multivariate normal distribution
         * at point `x`. There is no closed-form expression of this function. Hence it is evaluated
         * numerically via Monte-Carlo approximation.
         */

        if (!this->BelongsToSupport(x))
            throw itk::ExceptionObject(__FILE__, __LINE__, "The CDF is not defined outside the support.", ITK_LOCATION);

        unsigned int numParameters = x.size();

        if (m_MeanParameter.size() != numParameters)
            throw itk::ExceptionObject(__FILE__, __LINE__, "The input argument does not belong to the support as the one on which the normal distribution is defined.", ITK_LOCATION);

        const unsigned int numberOfMonteCarloSamples = 10000;
        SampleType normalSample(numberOfMonteCarloSamples);
        GeneratorType generator;
        this->Random(normalSample, generator);

        double sumValue = 0.0;
        for (unsigned int i = 0; i < numberOfMonteCarloSamples; ++i)
        {
            double sampleValue = 1.0;
            for (unsigned int j = 0; j < numParameters; ++j)
            {
                if (normalSample[i][j] > x[j])
                {
                    sampleValue = 0.0;
                    break;
                }
            }

            sumValue += sampleValue;
        }

        return sumValue / static_cast<double>(numberOfMonteCarloSamples);
    }

    void MultivariateNormalDistribution::Fit(const SampleType &sample, const std::string &method)
    {
        unsigned int numObservations = sample.size();
        unsigned int numParameters = m_MeanParameter.size();

        ValueType meanValue(numParameters);
        for (unsigned int j = 0; j < numParameters; ++j)
        {
            meanValue[j] = 0.0;
            for (unsigned int i = 0; i < numObservations; ++i)
                meanValue[j] += sample[i][j];
            meanValue[j] /= static_cast<double>(numObservations);
        }

        this->SetMeanParameter(meanValue);

        MatrixType covarianceMatrix(numParameters, numParameters);
        for (unsigned int j = 0; j < numParameters; ++j)
        {
            for (unsigned int k = 0; k < numParameters; ++k)
            {
                double matrixEntry = 0.0;
                for (unsigned int i = 0; i < numObservations; ++i)
                    matrixEntry += (sample[i][j] - meanValue[j]) * (sample[i][k] - meanValue[k]);
                matrixEntry /= static_cast<double>(numObservations);

                covarianceMatrix.put(j, k, matrixEntry);
                if (j != k)
                    covarianceMatrix.put(k, j, matrixEntry);
            }
        }

        this->SetCovarianceMatrixParameter(covarianceMatrix);
    }

    void MultivariateNormalDistribution::Random(SampleType &sample, GeneratorType &generator)
    {
        NormalDistributionType normDistr(0.0, 1.0);
        unsigned int numObservations = sample.size();
        unsigned int numParameters = m_MeanParameter.size();
        ValueType tmpValue(numParameters);

        for (unsigned int i = 0; i < numObservations; ++i)
        {
            for (unsigned int j = 0; j < numParameters; ++j)
                tmpValue[j] = normDistr(generator);

            sample[i].resize(numParameters);
            for (unsigned int j = 0; j < numParameters; ++j)
            {
                sample[i][j] = m_MeanParameter[j];
                for (unsigned int k = 0; k < numParameters; ++k)
                    sample[i][j] += m_SqrtCovarianceMatrixParameter.get(j, k) * tmpValue[k];
            }
        }
    }

    double MultivariateNormalDistribution::GetDistance(Self *otherDistribution)
    {
        /**
         * \fn double MultivariateNormalDistribution::GetDistance(MultivariateNormalDistribution *otherDistribution)
         *
         * \author Aymeric Stamm
         * \date November 2023
         *
         * \param otherDistribution A pointer specifying another object of class
         * `MultivariateNormalDistribution`.
         *
         * \return A numeric value storing the symmetric Kullback-Leibler divergence to the input normal
         * distribution. The implementation follows the formula in https://statproofbook.github.io/P/mvn-kl.
         */

        MultivariateNormalDistribution *normDistr = dynamic_cast<MultivariateNormalDistribution *>(otherDistribution);
        ValueType otherMeanParam = normDistr->GetMeanParameter();
        MatrixType otherCovParam = normDistr->GetCovarianceMatrixParameter();
        MatrixType otherPrecParam = normDistr->GetPrecisionMatrixParameter();

        unsigned int numParams = m_MeanParameter.size();
        if (otherMeanParam.size() != numParams)
            throw itk::ExceptionObject(__FILE__, __LINE__, "The two compared distributions should be on the same support.", ITK_LOCATION);

        double quadForm = 0.0;
        double firstTraceValue = 0.0, secondTraceValue = 0.0;
        for (unsigned int j = 0; j < numParams; ++j)
        {
            for (unsigned int k = 0; k < numParams; ++k)
            {
                quadForm += (otherMeanParam[j] - m_MeanParameter[j]) * (otherMeanParam[k] - m_MeanParameter[k]) * (m_PrecisionMatrixParameter.get(j, k) + otherPrecParam.get(j, k));
                firstTraceValue += otherPrecParam.get(j, k) * m_CovarianceMatrixParameter.get(k, j);
                secondTraceValue += m_PrecisionMatrixParameter.get(j, k) * otherCovParam.get(k, j);
            }
        }

        return quadForm / 2.0 + firstTraceValue + secondTraceValue - static_cast<double>(numParams);
    }
} // end of namespace anima
