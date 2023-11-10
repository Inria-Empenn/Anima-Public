#include "animaRealUniformDistribution.h"

#include <itkMacro.h>

namespace anima
{

    bool RealUniformDistribution::BelongsToSupport(const ValueType &x)
    {
        return (x >= m_LowerBoundParameter) && (x <= m_UpperBoundParameter);
    }

    double RealUniformDistribution::GetDensity(const ValueType &x)
    {
        if (!this->BelongsToSupport(x))
            return 0.0;
        return std::exp(this->GetLogDensity(x));
    }

    double RealUniformDistribution::GetLogDensity(const ValueType &x)
    {
        if (!this->BelongsToSupport(x))
            throw itk::ExceptionObject(__FILE__, __LINE__, "The log-density of the uniform distribution is not defined outside its support.", ITK_LOCATION);

        return -std::log(m_UpperBoundParameter - m_LowerBoundParameter);
    }

    double RealUniformDistribution::GetCumulative(const ValueType &x)
    {
        /**
         * \fn double RealUniformDistribution::GetCumulative(const double &x)
         *
         * \author Aymeric Stamm
         * \date November 2023
         *
         * \param x A numeric value specifying a point on the support of the real uniform
         * distribution.
         *
         * \return A numeric value storing the value of the cumulative distribution function
         * of the real uniform distribution at point `x`.
         */

        if (x < m_LowerBoundParameter)
            return 0.0;

        if (x > m_UpperBoundParameter)
            return 1.0;

        return (x - m_LowerBoundParameter) / (m_UpperBoundParameter - m_LowerBoundParameter);
    }

    void RealUniformDistribution::Fit(const SampleType &sample, const std::string &method)
    {
        this->SetLowerBoundParameter(*std::min_element(sample.begin(), sample.end()));
        this->SetUpperBoundParameter(*std::max_element(sample.begin(), sample.end()));
    }

    void RealUniformDistribution::Random(SampleType &sample, GeneratorType &generator)
    {
        UniformDistributionType unifDistr(m_LowerBoundParameter, m_UpperBoundParameter);
        unsigned int nSamples = sample.size();
        for (unsigned int i = 0; i < nSamples; ++i)
            sample[i] = unifDistr(generator);
    }

    double RealUniformDistribution::GetDistance(Self *otherDistribution)
    {
        /**
         * \fn double RealUniformDistribution::GetDistance(Self *otherDistribution)
         *
         * \author Aymeric Stamm
         * \date November 2023
         *
         * \param otherDistribution A pointer specifying another object of class `RealUniformDistribution`.
         *
         * \return A numeric value storing the symmetric Kullback-Leibler divergence with the
         * input real uniform distribution. See: https://statproofbook.github.io/P/cuni-kl.
         */

        return 0.0;
    }

} // end of namespace anima
