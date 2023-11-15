#include "animaSphericalUniformDistribution.h"
#include <animaVectorOperations.h>

#include <itkMacro.h>

namespace anima
{

    bool SphericalUniformDistribution::BelongsToSupport(const ValueType &x)
    {
        return std::abs(x.GetNorm() - 1.0) < this->GetEpsilon();
    }

    double SphericalUniformDistribution::GetDensity(const ValueType &x)
    {
        if (!this->BelongsToSupport(x))
            return 0.0;
        return std::exp(this->GetLogDensity(x));
    }

    double SphericalUniformDistribution::GetLogDensity(const ValueType &x)
    {
        if (!this->BelongsToSupport(x))
            throw itk::ExceptionObject(__FILE__, __LINE__, "The log-density is not defined outside the support.", ITK_LOCATION);

        return -std::log(4.0 * M_PI);
    }

    double SphericalUniformDistribution::GetCumulative(const ValueType &x)
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

        return phiVal / (2.0 * M_PI) * (1.0 - std::cos(thetaVal)) / 2.0;
    }

    void SphericalUniformDistribution::Random(SampleType &sample, GeneratorType &generator)
    {
        UniformDistributionType unifDistr(0.0, 1.0);
        ValueType sphCoords;
        sphCoords[2] = 1.0;
        unsigned int nSamples = sample.size();
        for (unsigned int i = 0; i < nSamples; ++i)
        {
            sphCoords[0] = 2.0 * std::asin(std::sqrt(unifDistr(generator)));
            sphCoords[1] = 2.0 * M_PI * unifDistr(generator);
            anima::TransformSphericalToCartesianCoordinates(sphCoords, sample[i]);
        }
    }

    SphericalUniformDistribution::ValueType SphericalUniformDistribution::GetMean()
    {
        ValueType meanValue;
        meanValue.Fill(0.0);
        return meanValue;
    }

    double SphericalUniformDistribution::GetDistance(Self *otherDistribution)
    {
        /**
         * \fn double SphericalUniformDistribution::GetDistance(Self *otherDistribution)
         *
         * \author Aymeric Stamm
         * \date November 2023
         *
         * \param otherDistribution A pointer specifying another object of class `SphericalUniformDistribution`.
         *
         * \return A numeric value storing the symmetric Kullback-Leibler divergence with the
         * input spherical uniform distribution.
         */

        return 0.0;
    }

} // end of namespace anima
