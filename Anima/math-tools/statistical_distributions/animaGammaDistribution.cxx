#include "animaGammaDistribution.h"
#include <animaGammaFunctions.h>

#include <itkMacro.h>

#include <boost/math/special_functions/gamma.hpp>

namespace anima
{

    bool GammaDistribution::BelongsToSupport(const ValueType &x)
    {
        return x > this->GetEpsilon();
    }

    void GammaDistribution::SetShapeParameter(const double val)
    {
        if (val < this->GetEpsilon())
            throw itk::ExceptionObject(__FILE__, __LINE__, "The shape parameter of the Gamma distribution should be strictly positive.", ITK_LOCATION);
        m_ShapeParameter = val;
    }

    void GammaDistribution::SetScaleParameter(const double val)
    {
        if (val < this->GetEpsilon())
            throw itk::ExceptionObject(__FILE__, __LINE__, "The scale parameter of the Gamma distribution should be strictly positive.", ITK_LOCATION);
        m_ScaleParameter = val;
    }

    double GammaDistribution::GetDensity(const ValueType &x)
    {
        if (!this->BelongsToSupport(x))
            return 0.0;
        return std::exp(this->GetLogDensity(x));
    }

    double GammaDistribution::GetLogDensity(const ValueType &x)
    {
        if (!this->BelongsToSupport(x))
            throw itk::ExceptionObject(__FILE__, __LINE__, "The log-density of the Gamma distribution is not defined for negative or null arguments.", ITK_LOCATION);
        ValueType resValue = (m_ShapeParameter - 1.0) * std::log(x);
        resValue -= x / m_ScaleParameter;
        resValue -= std::lgamma(m_ShapeParameter);
        resValue -= m_ShapeParameter * std::log(m_ScaleParameter);
        return resValue;
    }

    double GammaDistribution::GetCumulative(const ValueType &x)
    {
        /**
         * \fn double GetCumulative(const double &x)
         *
         * \author Aymeric Stamm (2023).
         *
         * \param x A numeric value specifying a point on the support of the Gamma distribution.
         *
         * \return A numeric value storing the value of the cumulative distribution function of
         * the Gamma distribution at point `x`. See: https://statproofbook.github.io/P/gam-cdf.
         */

        return boost::math::gamma_p<double, double>(this->GetShapeParameter(), x / this->GetScaleParameter());
    }

    void GammaDistribution::Fit(const SampleType &sample, const std::string &method)
    {
        unsigned int dimValue = sample.size();
        double doubleDimValue = static_cast<double>(dimValue);
        double epsValue = this->GetEpsilon();
        double shapeParameter, scaleParameter;

        if (method == "mle")
        {
            // Code for estimating theta and kappa via MLE
            double meanValue = 0.0;
            double meanLogValue = 0.0;
            for (unsigned int i = 0; i < dimValue; ++i)
            {
                double tmpValue = sample[i];
                if (tmpValue < epsValue)
                    throw itk::ExceptionObject(__FILE__, __LINE__, "Negative or null values are not allowed in a sample drawn from the Gamma distribution.", ITK_LOCATION);
                meanValue += tmpValue;
                meanLogValue += std::log(tmpValue);
            }
            meanValue /= doubleDimValue;
            meanLogValue /= doubleDimValue;
            double logMeanValue = std::log(meanValue);

            double sValue = logMeanValue - meanLogValue;
            shapeParameter = (3.0 - sValue + std::sqrt((sValue - 3.0) * (sValue - 3.0) + 24.0 * sValue)) / (12.0 * sValue);

            if (shapeParameter < epsValue)
                throw itk::ExceptionObject(__FILE__, __LINE__, "The shape parameter of a Gamma distribution cannot be negative or null.", ITK_LOCATION);

            scaleParameter = meanValue / shapeParameter;
        }
        else if (method == "biased-closed-form" || method == "unbiased-closed-form")
        {
            // Code for estimating theta and kappa via CF biased and unbiased
            double sumValue = 0.0;
            double sumLogValue = 0.0;
            double sumLogXValue = 0.0;
            for (unsigned int i = 0; i < dimValue; ++i)
            {
                double tmpValue = sample[i];
                if (tmpValue < std::numeric_limits<double>::epsilon())
                    throw itk::ExceptionObject(__FILE__, __LINE__, "Negative or null values are not allowed in a sample drawn from the Gamma distribution.", ITK_LOCATION);
                sumValue += tmpValue;
                sumLogValue += std::log(tmpValue);
                sumLogXValue += (std::log(tmpValue) * tmpValue);
            }

            double denValue = doubleDimValue * sumLogXValue - sumLogValue * sumValue;
            shapeParameter = doubleDimValue * sumValue;
            scaleParameter = doubleDimValue * sumLogXValue - sumLogValue * sumValue;
            if (denValue > epsValue)
                shapeParameter /= denValue;
            scaleParameter /= (doubleDimValue * doubleDimValue);

            if (method == "unbiased-closed-form")
            {
                shapeParameter -= (3.0 * shapeParameter - 2.0 / 3.0 * shapeParameter / (1.0 + shapeParameter) - 4.0 / 5.0 * shapeParameter / ((1.0 + shapeParameter) * (1.0 + shapeParameter))) / doubleDimValue;
                scaleParameter *= doubleDimValue / (doubleDimValue - 1.0);
            }
        }
        else
            throw itk::ExceptionObject(__FILE__, __LINE__, "Unsupported estimation method for the Gamma distribution.", ITK_LOCATION);

        this->SetShapeParameter(shapeParameter);
        this->SetScaleParameter(scaleParameter);
    }

    void GammaDistribution::Random(SampleType &sample, GeneratorType &generator)
    {
        DistributionType distributionValue(m_ShapeParameter, m_ScaleParameter);
        unsigned int nSamples = sample.size();
        for (unsigned int i = 0; i < nSamples; ++i)
            sample[i] = distributionValue(generator);
    }

    double GammaDistribution::GetDistance(Self *otherDistribution)
    {
        /**
         * \fn double GetDistance(GammaDistribution *otherDistribution)
         *
         * \author Aymeric Stamm (2023).
         *
         * \param otherDistribution A pointer specifying another object of class `GammaDistribution`.
         *
         * \return A numeric value storing the symmetric Kullback-Leibler divergence between the
         * current Gamma distribution and the input Gamma distribution. The calculation is done as
         * described in McCrimmon (2018), Distance metrics for Gamma distributions, arXiv:1802.01041v1.
         * Also documented here: https://statproofbook.github.io/P/gam-kl.
         */

        double thisKappa = this->GetShapeParameter();
        double thisTheta = this->GetScaleParameter();

        GammaDistribution *gammaDistr = dynamic_cast<GammaDistribution *>(otherDistribution);
        double otherKappa = gammaDistr->GetShapeParameter();
        double otherTheta = gammaDistr->GetScaleParameter();

        double distanceValue = (thisKappa - otherKappa) * (anima::digamma(thisKappa) + std::log(thisTheta) - anima::digamma(otherKappa) - std::log(otherTheta));
        distanceValue += (thisKappa * thisTheta - otherKappa * otherTheta) * (thisTheta - otherTheta) / (thisTheta * otherTheta);
        return distanceValue;
    }

} // end of namespace anima
