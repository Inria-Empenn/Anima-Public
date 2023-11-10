#include "animaWatsonDistribution.h"
#include <animaErrorFunctions.h>
#include <animaKummerFunctions.h>
#include <animaMatrixOperations.h>

#include <itkMacro.h>
#include <itkSymmetricEigenAnalysis.h>

namespace anima
{

    bool WatsonDistribution::BelongsToSupport(const ValueType &x)
    {
        return std::abs(x.GetNorm() - 1.0) < this->GetEpsilon();
    }

    void WatsonDistribution::SetMeanAxis(const ValueType &val)
    {
        if (!this->BelongsToSupport(val))
            throw itk::ExceptionObject(__FILE__, __LINE__, "The mean axis parameter of the Watson distribution should be of unit norm.", ITK_LOCATION);
        m_MeanAxis = val;
    }

    void WatsonDistribution::SetConcentrationParameter(const double &val)
    {
        m_ConcentrationParameter = val;
    }

    double WatsonDistribution::GetDensity(const ValueType &x)
    {
        if (!this->BelongsToSupport(x))
            return 0.0;

        // Case 1: k = 0 - Uniform distribution on the 2-sphere
        if (std::abs(m_ConcentrationParameter) <= this->GetEpsilon())
            return 1.0 / (4.0 * M_PI);

        // Case 2: k > 0 - Anisotropic distribution on the 2-sphere
        if (m_ConcentrationParameter > this->GetEpsilon())
        {
            double kappaSqrt = std::sqrt(m_ConcentrationParameter);
            double c = anima::ComputeScalarProduct(x, m_MeanAxis);
            double inExp = m_ConcentrationParameter * (c * c - 1.0);
            return kappaSqrt * std::exp(inExp) / (4.0 * M_PI * anima::EvaluateDawsonFunctionNR(kappaSqrt));
        }

        // Case 3: k < 0 - Gridle distribution on the 2-sphere
        double Ck = std::sqrt(-m_ConcentrationParameter / M_PI) / (2.0 * M_PI * std::erf(std::sqrt(-m_ConcentrationParameter)));
        double c = anima::ComputeScalarProduct(x, m_MeanAxis);
        double inExp = m_ConcentrationParameter * c * c;
        return Ck * std::exp(inExp);
    }

    double WatsonDistribution::GetLogDensity(const ValueType &x)
    {
        if (!this->BelongsToSupport(x))
            throw itk::ExceptionObject(__FILE__, __LINE__, "The log-density of the Watson distribution is not defined for arguments outside the 2-sphere.", ITK_LOCATION);

        return std::log(this->GetDensity(x));
    }

    double WatsonDistribution::ComputeConcentrationMLE(const double rValue, const double aValue, const double cValue, double &logLik)
    {
        double bBoundValue = (rValue * cValue - aValue) / (2.0 * rValue * (1.0 - rValue));
        double rootInValue = 1.0 + (4.0 * (cValue + 1.0) * rValue * (1.0 - rValue)) / (aValue * (cValue - aValue));
        bBoundValue *= (1.0 + std::sqrt(rootInValue));

        double lBoundValue = (rValue * cValue - aValue) / (rValue * (1.0 - rValue));
        lBoundValue *= (1.0 + (1.0 - rValue) / (cValue - aValue));

        double uBoundValue = (rValue * cValue - aValue) / (rValue * (1.0 - rValue));
        uBoundValue *= (1.0 + rValue / aValue);

        double concentrationParameter = 0.0;
        if (rValue > aValue / cValue)
            concentrationParameter = (bBoundValue + lBoundValue) / 2.0;
        if (rValue < aValue / cValue)
            concentrationParameter = (bBoundValue + uBoundValue) / 2.0;

        logLik = concentrationParameter * (rValue - 1.0) - std::log(anima::GetScaledKummerFunctionValue(concentrationParameter, aValue, cValue));

        return concentrationParameter;
    }

    void WatsonDistribution::Fit(const SampleType &sample, const std::string &method)
    {
        /**********************************************************************************************
         * \fn void Fit(std::vector<itk::Vector<double,3>>, std::mt19937 &generator)
         *
         * \brief	Closed-form approximations of the maximum likelihood estimators for the mean axis
         *          and the concentration parameter of the Watson distribution using the procedure
         *          described in Sra & Karp, The multivariate Watson distribution: Maximum-likelihood
         *          estimation and other aspects, Journal of Multivariate Analysis, 2013, Vol. 114,
         *          pp. 256-69.
         *
         * \author	Aymeric Stamm
         * \date	October 2023
         *
         * \param	sample A numeric matrix of shape n x 3 specifying a sample of size n drawn from the
         *          Watson distribution on the 2-sphere.
         * \param	method A string specifying the estimation method. Unused here.
         **********************************************************************************************/

        unsigned int numberOfObservations = sample.size();
        ValueType meanAxis;
        meanAxis.Fill(0.0);
        meanAxis[2] = 1.0;
        double concentrationParameter = 0.0;

        using MatrixType = itk::Matrix<double, 3, 3>;
        MatrixType scatterMatrix;
        scatterMatrix.Fill(0.0);
        for (unsigned int i = 0; i < numberOfObservations; ++i)
        {
            for (unsigned int j = 0; j < m_AmbientDimension; ++j)
            {
                for (unsigned int k = j; k < m_AmbientDimension; ++k)
                {
                    double tmpValue = sample[i][j] * sample[i][k];
                    scatterMatrix(j, k) += tmpValue;
                    if (j != k)
                        scatterMatrix(k, j) += tmpValue;
                }
            }
        }

        scatterMatrix /= static_cast<double>(numberOfObservations);

        // Eigen decomposition of S
        itk::SymmetricEigenAnalysis<MatrixType, ValueType, MatrixType> eigenDecomp(m_AmbientDimension);
        ValueType eigenVals;
        MatrixType eigenVecs;
        eigenDecomp.ComputeEigenValuesAndVectors(scatterMatrix, eigenVals, eigenVecs);

        double aValue = 0.5;
        double cValue = static_cast<double>(m_AmbientDimension) / 2.0;

        // Compute estimate of mean axis assuming estimated concentration parameter is positive
        double bipolarRValue = 0.0;
        for (unsigned int i = 0; i < m_AmbientDimension; ++i)
        {
            for (unsigned int j = 0; j < m_AmbientDimension; ++j)
                bipolarRValue += eigenVecs(2, i) * scatterMatrix(i, j) * eigenVecs(2, j);
        }

        double bipolarLogLik;
        double bipolarKappa = ComputeConcentrationMLE(bipolarRValue, aValue, cValue, bipolarLogLik);
        bool bipolarValid = bipolarKappa > 0;

        // Compute estimate of mean axis assuming estimated concentration parameter is negative
        double gridleRValue = 0.0;
        for (unsigned int i = 0; i < m_AmbientDimension; ++i)
        {
            for (unsigned int j = 0; j < m_AmbientDimension; ++j)
                gridleRValue += eigenVecs(0, i) * scatterMatrix(i, j) * eigenVecs(0, j);
        }

        double gridleLogLik;
        double gridleKappa = ComputeConcentrationMLE(gridleRValue, aValue, cValue, gridleLogLik);
        bool gridleValid = gridleKappa < 0;

        if (bipolarValid && gridleValid)
        {
            if (gridleLogLik > bipolarLogLik)
            {
                for (unsigned int i = 0; i < m_AmbientDimension; ++i)
                    meanAxis[i] = eigenVecs(0, i);
                concentrationParameter = gridleKappa;
                m_RValue = gridleRValue;
            }
            else
            {
                for (unsigned int i = 0; i < m_AmbientDimension; ++i)
                    meanAxis[i] = eigenVecs(2, i);
                concentrationParameter = bipolarKappa;
                m_RValue = bipolarRValue;
            }
        }
        else if (bipolarValid)
        {
            for (unsigned int i = 0; i < m_AmbientDimension; ++i)
                meanAxis[i] = eigenVecs(2, i);
            concentrationParameter = bipolarKappa;
            m_RValue = bipolarRValue;
        }
        else if (gridleValid)
        {
            for (unsigned int i = 0; i < m_AmbientDimension; ++i)
                meanAxis[i] = eigenVecs(0, i);
            concentrationParameter = bipolarKappa;
            m_RValue = gridleRValue;
        }
        else
        {
            concentrationParameter = 0.0;
            m_RValue = gridleRValue;
        }

        meanAxis.Normalize();

        this->SetMeanAxis(meanAxis);
        this->SetConcentrationParameter(concentrationParameter);
    }

    void WatsonDistribution::Random(SampleType &sample, GeneratorType &generator)
    {
        /**********************************************************************************************
         * \fn void Random(std::vector<itk::Vector<double,3>>, std::mt19937 &generator)
         *
         * \brief	Sample from the Watson distribution using the procedure described in Fisher et al.,
         *          Statistical Analysis of Spherical Data, Cambridge University Press, 1993, pp. 59.
         *
         * \author	Aymeric Stamm
         * \date	October 2013
         *
         * \param	sample    A numeric matrix of shape n x 3 storing a sample of size n drawn from the
         *                    Watson distribution on the 2-sphere.
         * \param	generator A pseudo-random number generator.
         **********************************************************************************************/

        ValueType tmpVec, resVec;
        tmpVec.Fill(0.0);
        tmpVec[2] = 1.0;

        // Compute rotation matrix to bring [0,0,1] on meanAxis
        itk::Matrix<double, 3, 3> rotationMatrix = anima::GetRotationMatrixFromVectors(tmpVec, m_MeanAxis);

        // Now resuming onto sampling around direction [0,0,1]
        double U, V, S;
        UniformDistributionType distributionValue(0.0, 1.0);
        unsigned int nSamples = sample.size();

        for (unsigned int i = 0; i < nSamples; ++i)
        {
            if (m_ConcentrationParameter > std::sqrt(std::numeric_limits<double>::epsilon())) // Bipolar distribution
            {
                U = distributionValue(generator);
                S = 1.0 + std::log(U + (1.0 - U) * std::exp(-m_ConcentrationParameter)) / m_ConcentrationParameter;

                V = distributionValue(generator);

                if (V > 1.0e-6)
                {
                    while (std::log(V) > m_ConcentrationParameter * S * (S - 1.0))
                    {
                        U = distributionValue(generator);
                        S = 1.0 + std::log(U + (1.0 - U) * std::exp(-m_ConcentrationParameter)) / m_ConcentrationParameter;

                        V = distributionValue(generator);

                        if (V < 1.0e-6)
                            break;
                    }
                }
            }
            else if (m_ConcentrationParameter < -std::sqrt(std::numeric_limits<double>::epsilon())) // Gridle distribution
            {
                double C1 = std::sqrt(std::abs(m_ConcentrationParameter));
                double C2 = std::atan(C1);
                U = distributionValue(generator);
                V = distributionValue(generator);
                S = (1.0 / C1) * std::tan(C2 * U);

                double T = m_ConcentrationParameter * S * S;
                while (V > (1.0 - T) * std::exp(T))
                {
                    U = distributionValue(generator);
                    V = distributionValue(generator);
                    S = (1.0 / C1) * std::tan(C2 * U);
                    T = m_ConcentrationParameter * S * S;
                }
            }
            else // Uniform distribution
                S = std::cos(M_PI * distributionValue(generator));

            double phi = 2.0 * M_PI * distributionValue(generator);

            tmpVec[0] = std::sqrt(1.0 - S * S) * std::cos(phi);
            tmpVec[1] = std::sqrt(1.0 - S * S) * std::sin(phi);
            tmpVec[2] = S;

            // Rotate to bring everthing back around meanAxis
            for (unsigned int j = 0; j < m_AmbientDimension; ++j)
            {
                resVec[j] = 0.0;
                for (unsigned int k = 0; k < m_AmbientDimension; ++k)
                    resVec[j] += rotationMatrix(j, k) * tmpVec[k];
            }

            double resNorm = resVec.Normalize();

            if (std::abs(resNorm - 1.0) > this->GetEpsilon())
            {
                std::cout << "Sampled direction norm: " << resNorm << std::endl;
                std::cout << "Mean direction: " << m_MeanAxis << std::endl;
                std::cout << "Concentration parameter: " << m_ConcentrationParameter << std::endl;
                throw itk::ExceptionObject(__FILE__, __LINE__, "The Watson sampler should generate points on the 2-sphere.", ITK_LOCATION);
            }

            sample[i] = resVec;
        }
    }

    void WatsonDistribution::GetStandardWatsonSHCoefficients(std::vector<double> &coefficients,
                                                             std::vector<double> &derivatives)
    {
        // Computes the first 7 non-zero SH coefficients of the standard Watson PDF (multiplied by 4 M_PI).
        const unsigned int nbCoefs = 7;
        coefficients.resize(nbCoefs);
        derivatives.resize(nbCoefs);

        double sqrtPi = std::sqrt(M_PI);
        double k = m_ConcentrationParameter;
        double dawsonValue = anima::EvaluateDawsonIntegral(std::sqrt(k), true);
        double k2 = k * k;
        double k3 = k2 * k;
        double k4 = k3 * k;
        double k5 = k4 * k;
        double k6 = k5 * k;
        double k7 = k6 * k;

        coefficients[0] = 2.0 * sqrtPi;
        derivatives[0] = 0.0;

        if (k < 1.0e-4) // Handles small k values by Taylor series expansion
        {
            coefficients[1] = k / 3.0;
            coefficients[2] = k2 / 35.0;
            coefficients[3] = k3 / 693.0;
            coefficients[4] = k4 / 19305.0;
            coefficients[5] = k5 / 692835.0;
            coefficients[6] = k6 / 30421755.0;
            derivatives[1] = 1.0 / 3.0;
            derivatives[2] = 2.0 * k / 35.0;
            derivatives[3] = k2 / 231.0;
            derivatives[4] = 4.0 * k3 / 19305.0;
            derivatives[5] = k4 / 138567.0;
            derivatives[6] = 6.0 * k5 / 30421755.0;

            for (unsigned int i = 1; i < nbCoefs; ++i)
            {
                double tmpVal = std::pow(2.0, i + 1.0) * sqrtPi / std::sqrt(4.0 * i + 1.0);
                coefficients[i] *= tmpVal;
                derivatives[i] *= tmpVal;
            }

            return;
        }

        coefficients[1] = (3.0 - (3.0 + 2.0 * k) * dawsonValue) / (2.0 * k);
        coefficients[2] = (5.0 * (-21.0 + 2.0 * k) + 3.0 * (35.0 + 4.0 * k * (5.0 + k)) * dawsonValue) / (16.0 * k2);
        coefficients[3] = (21.0 * (165.0 + 4.0 * (-5.0 + k) * k) - 5.0 * (693.0 + 378.0 * k + 84.0 * k2 + 8.0 * k3) * dawsonValue) / (64.0 * k3);
        coefficients[4] = (3.0 * (-225225.0 + 2.0 * k * (15015.0 + 2.0 * k * (-1925.0 + 62.0 * k))) + 35.0 * (19305.0 + 8.0 * k * (1287.0 + k * (297.0 + 2.0 * k * (18.0 + k)))) * dawsonValue) / (1024.0 * k4);
        coefficients[5] = (11.0 * (3968055.0 + 8.0 * k * (-69615.0 + 2.0 * k * (9828.0 + k * (-468.0 + 29.0 * k)))) - 63.0 * (692835.0 + 2.0 * k * (182325.0 + 4.0 * k * (10725.0 + 2.0 * k * (715.0 + k * (55.0 + 2.0 * k))))) * dawsonValue) / (4096.0 * k5);
        coefficients[6] = (13.0 * (-540571185.0 + 2.0 * k * (39171825.0 + 4.0 * k * (-2909907.0 + 2.0 * k * (82467.0 + k * (-7469.0 + 122.0 * k))))) + 231.0 * (30421755.0 + 4.0 * k * (3968055.0 + k * (944775.0 + 4.0 * k * (33150.0 + k * (2925.0 + 4.0 * k * (39.0 + k)))))) * dawsonValue) / (32768.0 * k6);

        derivatives[1] = 3.0 * (-1.0 + (-1.0 + 2.0 * k) * dawsonValue + 2.0 * dawsonValue * dawsonValue) / (4.0 * k2);
        derivatives[2] = 5.0 * ((21.0 - 2.0 * k) + (63.0 + 4.0 * (-11.0 + k) * k) * dawsonValue - 12.0 * (7.0 + 2.0 * k) * dawsonValue * dawsonValue) / (32.0 * k3);
        derivatives[3] = 21.0 * ((-165.0 - 4.0 * (-5.0 + k) * k) + (-5.0 + 2.0 * k) * (165.0 + 4.0 * (-3.0 + k) * k) * dawsonValue + 10.0 * (99.0 + 4.0 * k * (9.0 + k)) * dawsonValue * dawsonValue) / (128.0 * k4);
        derivatives[4] = 3.0 * ((225225.0 - 2.0 * k * (15015.0 + 2.0 * k * (-1925.0 + 62.0 * k))) + (1576575.0 + 8.0 * k * (-75075.0 + k * (10395.0 + 2.0 * k * (-978.0 + 31.0 * k)))) * dawsonValue - 840.0 * (2145.0 + 2.0 * k * (429.0 + 66.0 * k + 4.0 * k2)) * dawsonValue * dawsonValue) / (2048.0 * k5);
        derivatives[5] = 11.0 * ((-3968055.0 - 8.0 * k * (-69615.0 + 2.0 * k * (9828.0 + k * (-468.0 + 29.0 * k)))) + (-35712495.0 + 2.0 * k * (5917275.0 + 8.0 * k * (-118755.0 + k * (21060.0 + k * (-965.0 + 58.0 * k))))) * dawsonValue + 630.0 * (62985.0 + 8.0 * k * (3315.0 + k * (585.0 + 2.0 * k * (26.0 + k)))) * dawsonValue * dawsonValue) / (8192.0 * k6);
        derivatives[6] = 13.0 * ((540571185.0 - 2.0 * k * (39171825.0 + 4.0 * k * (-2909907.0 + 2.0 * k * (82467.0 + k * (-7469.0 + 122.0 * k))))) + (5946283035.0 + 4.0 * k * (-446558805.0 + k * (79910523.0 + 4.0 * k * (-3322242.0 + k * (187341.0 + 4.0 * k * (-3765.0 + 61.0 * k)))))) * dawsonValue - 2772.0 * (2340135.0 + 2.0 * k * (508725.0 + 4.0 * k * (24225.0 + 2.0 * k * (1275.0 + k * (75.0 + 2.0 * k))))) * dawsonValue * dawsonValue) / (65536.0 * k7);

        for (unsigned int i = 1; i < nbCoefs; ++i)
        {
            double sqrtVal = std::sqrt(1.0 + 4.0 * i);
            coefficients[i] *= (sqrtPi * sqrtVal / dawsonValue);
            derivatives[i] *= (sqrtPi * sqrtVal / (dawsonValue * dawsonValue));
        }
    }

    double WatsonDistribution::GetDistance(Self *otherDistribution)
    {
        const unsigned int numberOfMonteCarloSamples = 10000;
        SampleType thisWatsonSample(numberOfMonteCarloSamples), otherWatsonSample(numberOfMonteCarloSamples);
        GeneratorType generator;

        this->Random(thisWatsonSample, generator);
        WatsonDistribution *watsonDistr = dynamic_cast<WatsonDistribution *>(otherDistribution);
        watsonDistr->Random(otherWatsonSample, generator);

        double thisKLValue = 0.0, otherKLValue = 0.0;
        for (unsigned int i = 0; i < numberOfMonteCarloSamples; ++i)
        {
            thisKLValue += this->GetLogDensity(thisWatsonSample[i]);
            thisKLValue -= watsonDistr->GetLogDensity(thisWatsonSample[i]);
            otherKLValue += watsonDistr->GetLogDensity(otherWatsonSample[i]);
            otherKLValue -= this->GetLogDensity(otherWatsonSample[i]);
        }

        thisKLValue /= static_cast<double>(numberOfMonteCarloSamples);
        otherKLValue /= static_cast<double>(numberOfMonteCarloSamples);

        return thisKLValue + otherKLValue;
    }

} // end of namespace anima
