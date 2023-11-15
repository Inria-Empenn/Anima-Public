#include "animaVonMisesFisherDistribution.h"
#include <animaBesselFunctions.h>
#include <animaMatrixOperations.h>
#include <animaVectorOperations.h>

#include <itkMacro.h>

namespace anima
{

    bool VonMisesFisherDistribution::BelongsToSupport(const ValueType &x)
    {
        return std::abs(x.GetNorm() - 1.0) < this->GetEpsilon();
    }

    void VonMisesFisherDistribution::SetMeanDirection(const ValueType &val)
    {
        if (!this->BelongsToSupport(val))
            throw itk::ExceptionObject(__FILE__, __LINE__, "The mean direction parameter of the von Mises Fisher distribution should be of unit norm.", ITK_LOCATION);
        m_MeanDirection = val;
    }

    void VonMisesFisherDistribution::SetConcentrationParameter(const double &val)
    {
        if (val < this->GetEpsilon())
            throw itk::ExceptionObject(__FILE__, __LINE__, "The concentration parameter of the von Mises Fisher distribution should be positive.", ITK_LOCATION);
        m_ConcentrationParameter = val;
        m_BesselRatio = anima::bessel_ratio_i_lower_bound(val, static_cast<double>(m_AmbientDimension) / 2.0);
    }

    double VonMisesFisherDistribution::GetDensity(const ValueType &x)
    {
        if (!this->BelongsToSupport(x))
            return 0.0;

        double kappa = this->GetConcentrationParameter();
        ValueType meanDirection = this->GetMeanDirection();

        if (kappa < std::sqrt(this->GetEpsilon()))
            return std::exp(kappa * anima::ComputeScalarProduct(meanDirection, x)) / (4.0 * M_PI);

        double tmpVal = kappa * (anima::ComputeScalarProduct(meanDirection, x) - 1.0);
        double resVal = std::exp(tmpVal);
        resVal *= kappa;
        tmpVal = 1.0 - std::exp(-2.0 * kappa);
        resVal /= (2.0 * M_PI * tmpVal);

        return resVal;
    }

    double VonMisesFisherDistribution::GetLogDensity(const ValueType &x)
    {
        if (!this->BelongsToSupport(x))
            throw itk::ExceptionObject(__FILE__, __LINE__, "The log-density of the von Mises Fisher distribution is not defined for arguments outside the 2-sphere.", ITK_LOCATION);

        return std::log(this->GetDensity(x));
    }

    double VonMisesFisherDistribution::GetCumulative(const ValueType &x)
    {
        if (!this->BelongsToSupport(x))
            throw itk::ExceptionObject(__FILE__, __LINE__, "The CDF is not defined outside the support.", ITK_LOCATION);

        ValueType sphCoords;
        anima::TransformCartesianToSphericalCoordinates(x, sphCoords);
        double thetaVal = sphCoords[0];
        while (thetaVal > M_PI)
            thetaVal -= (2.0 * M_PI);
        while (thetaVal < 0)
            thetaVal += (2.0 * M_PI);
        double phiVal = sphCoords[1];
        while (phiVal > 2.0 * M_PI)
            phiVal -= (2.0 * M_PI);
        while (phiVal < 0)
            phiVal += 2.0 * M_PI;

        double phiCumul = phiVal / (2.0 * M_PI);
        double thetaCumul = (1.0 - std::exp(-m_ConcentrationParameter * (1.0 - std::cos(thetaVal))));
        thetaCumul /= (1.0 - std::exp(-2.0 * m_ConcentrationParameter));

        return phiCumul * thetaCumul;
    }

    void VonMisesFisherDistribution::Fit(const SampleType &sample, const std::string &method)
    {
        /**********************************************************************************************
         * \fn      void VonMisesFisherDistribution::Fit(std::vector<itk::Vector<double,3>>,
         *                                               std::mt19937 &generator)
         *
         * \brief	Closed-form approximations of the maximum likelihood estimators for the mean
         *          direction and concentration parameter of the von Mises Fisher distribution using
         *          the procedure described in Sra, A short note on parameter approximation for von
         *          Mises-Fisher distributions: and a fast implementation of Is(x), Computational
         *          Statistics, 2011.
         *
         * \author	Aymeric Stamm
         * \date	October 2023
         *
         * \param	sample A numeric matrix of shape n x 3 specifying a sample of size n drawn from the
         *          von Mises Fisher distribution on the 2-sphere.
         * \param	method A string specifying the estimation method. Unused here.
         **********************************************************************************************/

        unsigned int numberOfObservations = sample.size();
        double ambientDimension = static_cast<double>(m_AmbientDimension);

        // Eq. (2)
        ValueType meanDirection;
        meanDirection.Fill(0.0);
        for (unsigned int i = 0; i < numberOfObservations; ++i)
            for (unsigned int j = 0; j < m_AmbientDimension; ++j)
                meanDirection[j] += sample[i][j];
        double normValue = meanDirection.Normalize();
        double resultantValue = normValue / static_cast<double>(numberOfObservations);

        if (std::abs(resultantValue - 1.0) < this->GetEpsilon())
        {
            // All observations are perfectly aligned
            // Kappa should be close to 0.
            this->SetMeanDirection(meanDirection);
            this->SetConcentrationParameter(this->GetEpsilon());
            return;
        }

        // Eq. (4)
        double concentrationParameter = resultantValue * (ambientDimension - resultantValue * resultantValue);
        concentrationParameter /= (1.0 - resultantValue * resultantValue);

        // Eq. (6)
        bool continueLoop = true;
        while (continueLoop)
        {
            double oldConcentrationParameter = concentrationParameter;
            double besselRatio = anima::bessel_ratio_i_lower_bound(concentrationParameter, ambientDimension / 2.0);
            double tmpValue = besselRatio - resultantValue;
            tmpValue /= (1.0 - besselRatio * besselRatio - (ambientDimension - 1.0) / concentrationParameter * besselRatio);
            concentrationParameter -= tmpValue;
            if (std::abs(concentrationParameter - oldConcentrationParameter) < this->GetEpsilon())
                continueLoop = false;
        }

        this->SetMeanDirection(meanDirection);
        this->SetConcentrationParameter(concentrationParameter);
    }

    void VonMisesFisherDistribution::Random(SampleType &sample, GeneratorType &generator)
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

        unsigned int sampleSize = sample.size();
        ValueType sampleValue;

        if (m_ConcentrationParameter > 700.0)
        {
            for (unsigned int i = 0; i < sampleSize; ++i)
            {
                this->SampleFromVMFDistributionNumericallyStable(sampleValue, generator);
                sample[i] = sampleValue;
            }
        }
        else
        {
            for (unsigned int i = 0; i < sampleSize; ++i)
            {
                this->SampleFromVMFDistribution(sampleValue, generator);
                sample[i] = sampleValue;
            }
        }
    }

    double VonMisesFisherDistribution::GetDistance(Self *otherDistribution)
    {
        /**
         * \fn double VonMisesFisherDistribution::GetDistance(VonMisesFisherDistribution *otherDistribution)
         *
         * \author Aymeric Stamm
         * \date November 2023
         *
         * \param otherDistribution An object of class `VonMisesFisherDistribution`.
         *
         * \return A numeric value storing the symmetric Kullback-Leibler divergence to the input von Mises
         * Fisher distribution. This is achieved following Kitagawa & Rowley (2022), von Mises-Fisher
         * distributions and their statistical divergence, arXiv:2202.05192v1
         * (https://arxiv.org/pdf/2202.05192v1.pdf).
         */

        VonMisesFisherDistribution *vmfDistr = dynamic_cast<VonMisesFisherDistribution *>(otherDistribution);
        ValueType otherMeanDirection = vmfDistr->GetMeanDirection();

        if (otherMeanDirection.GetVectorDimension() != m_AmbientDimension)
            throw itk::ExceptionObject(__FILE__, __LINE__, "The two compared distributions should be on the same sphere.", ITK_LOCATION);

        double otherConcentrationParameter = vmfDistr->GetConcentrationParameter();
        double otherBesselRatio = vmfDistr->GetBesselRatio();

        double thisToOtherDist = 0.0;
        for (unsigned int i = 0; i < m_AmbientDimension; ++i)
            thisToOtherDist += (m_ConcentrationParameter * m_MeanDirection[i] - otherConcentrationParameter * otherMeanDirection[i]) * m_MeanDirection[i];
        thisToOtherDist *= m_BesselRatio;

        double otherToThisDist = 0.0;
        for (unsigned int i = 0; i < m_AmbientDimension; ++i)
            otherToThisDist += (otherConcentrationParameter * otherMeanDirection[i] - m_ConcentrationParameter * m_MeanDirection[i]) * otherMeanDirection[i];
        otherToThisDist *= otherBesselRatio;

        return thisToOtherDist + otherToThisDist;
    }

    vnl_matrix<double> VonMisesFisherDistribution::GetCovarianceMatrix()
    {
        /**
         * \fn vnl_matrix<double> VonMisesFisherDistribution::GetCovarianceMatrix()
         *
         * \author Aymeric Stamm
         * \date November 2023
         *
         * \return A numeric matrix of size `m_AmbientDimension x m_AmbientDimension` storing the
         * covariance matrix of the von Mises Fisher distribution. This is achieved following
         * Kitagawa & Rowley (2022), von Mises-Fisher distributions and their statistical divergence,
         * arXiv:2202.05192v1 (https://arxiv.org/pdf/2202.05192v1.pdf).
         */

        vnl_matrix<double> covarianceMatrix(m_AmbientDimension, m_AmbientDimension);
        double diagConstant = m_BesselRatio / m_ConcentrationParameter;
        double offDiagConstant = 1.0 - static_cast<double>(m_AmbientDimension) * m_BesselRatio / m_ConcentrationParameter - m_BesselRatio * m_BesselRatio;

        for (unsigned int i = 0; i < m_AmbientDimension; ++i)
        {
            covarianceMatrix(i, i) = diagConstant + offDiagConstant * m_MeanDirection[i] * m_MeanDirection[i];

            for (unsigned int j = i + 1; j < m_AmbientDimension; ++j)
            {
                double tmpValue = offDiagConstant * m_MeanDirection[i] * m_MeanDirection[j];
                covarianceMatrix(i, j) = tmpValue;
                covarianceMatrix(j, i) = tmpValue;
            }
        }

        return covarianceMatrix;
    }

    void VonMisesFisherDistribution::SampleFromVMFDistribution(ValueType &resVec, GeneratorType &generator)
    {
        /**
         * \fn void VonMisesFisherDistribution::SampleFromVMFDistribution(&resVec, &generator)
         *
         * \brief Sample from the VMF distribution following Ulrich, G. (1984). Computer generation of
         * distributions on the m‐sphere. Journal of the Royal Statistical Society: Series C (Applied
         * Statistics), 33(2), 158-163.
         *
         * \author Aymeric Stamm
         * \date October 2013
         *
         * \param resVec An object of type itk::Vector<double,3> that will store a value sampled from
         * the VMF distribution.
         * \param generator An object of type std::mt19937 specifying which random number generator to
         * use.
         */

        ValueType tmpVec;
        for (unsigned int i = 0; i < m_AmbientDimension; ++i)
            tmpVec[i] = 0;
        tmpVec[2] = 1.0;

        // Compute rotation matrix to bring [0,0,1] on meanDirection
        itk::Matrix<double, 3, 3> rotationMatrix = anima::GetRotationMatrixFromVectors(tmpVec, m_MeanDirection);

        // Now resuming onto sampling around direction [0,0,1]
        UniformDistributionType unifDistr(0.0, 1.0);
        BetaDistributionType betaDistr(1.0, 1.0);

        double tmpVal = std::sqrt(m_ConcentrationParameter * m_ConcentrationParameter + 1.0);
        double b = (-2.0 * m_ConcentrationParameter + 2.0 * tmpVal) / 2.0;
        double a = (1.0 + m_ConcentrationParameter + tmpVal) / 2.0;
        double d = 4.0 * a * b / (1.0 + b) - 2.0 * std::log(2.0);

        double T = 1.0;
        double U = std::exp(d);
        double W = 0;

        while (2.0 * std::log(T) - T + d < std::log(U))
        {
            double Z = boost::math::quantile(betaDistr, unifDistr(generator));
            U = unifDistr(generator);
            tmpVal = 1.0 - (1.0 - b) * Z;
            T = 2.0 * a * b / tmpVal;
            W = (1.0 - (1.0 + b) * Z) / tmpVal;
        }

        double theta = 2.0 * M_PI * unifDistr(generator);
        tmpVec[0] = std::sqrt(1.0 - W * W) * std::cos(theta);
        tmpVec[1] = std::sqrt(1.0 - W * W) * std::sin(theta);
        tmpVec[2] = W;

        // Rotate to bring everthing back around meanDirection
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
            std::cout << "Mean direction: " << m_MeanDirection << std::endl;
            std::cout << "Concentration parameter: " << m_ConcentrationParameter << std::endl;
            throw itk::ExceptionObject(__FILE__, __LINE__, "The VMF sampler should generate points on the 2-sphere.", ITK_LOCATION);
        }
    }

    void VonMisesFisherDistribution::SampleFromVMFDistributionNumericallyStable(ValueType &resVec, GeneratorType &generator)
    {
        /**
         * \fn void VonMisesFisherDistribution::SampleFromVMFDistributionNumericallyStable(&resVec, &generator)
         *
         * \brief Sample from the VMF distribution following Jakob, W. (2012). Numerically stable sampling
         * of the von Mises-Fisher distribution on Sˆ2 (and other tricks). Interactive Geometry Lab, ETH
         * Zürich, Tech. Rep, 6.
         *
         * \author Aymeric Stamm
         * \date October 2013
         *
         * \param resVec An object of type itk::Vector<double,3> that will store a value sampled from
         * the VMF distribution.
         * \param generator An object of type std::mt19937 specifying which random number generator to
         * use.
         */

        ValueType tmpVec;
        for (unsigned int i = 0; i < m_AmbientDimension; ++i)
            tmpVec[i] = 0;
        tmpVec[2] = 1.0;

        // Compute rotation matrix to bring [0,0,1] on meanDirection
        itk::Matrix<double, 3, 3> rotationMatrix = anima::GetRotationMatrixFromVectors(tmpVec, m_MeanDirection);

        // Now resuming onto sampling around direction [0,0,1]
        UniformDistributionType unifDistr(0.0, 1.0);

        double xi = unifDistr(generator);
        double W = 1.0 + (std::log(xi) + std::log(1.0 - (xi - 1.0) * exp(-2.0 * m_ConcentrationParameter) / xi)) / m_ConcentrationParameter;
        double theta = 2.0 * M_PI * unifDistr(generator);

        tmpVec[0] = std::sqrt(1.0 - W * W) * std::cos(theta);
        tmpVec[1] = std::sqrt(1.0 - W * W) * std::sin(theta);
        tmpVec[2] = W;

        // Rotate to bring everthing back around meanDirection
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
            std::cout << "Mean direction: " << m_MeanDirection << std::endl;
            std::cout << "Concentration parameter: " << m_ConcentrationParameter << std::endl;
            throw itk::ExceptionObject(__FILE__, __LINE__, "The VMF sampler should generate points on the 2-sphere.", ITK_LOCATION);
        }
    }

} // end of namespace anima
