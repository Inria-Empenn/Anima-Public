#pragma once

#include <animaBaseDistribution.h>

namespace anima
{
    class ANIMASTATISTICALDISTRIBUTIONS_EXPORT RealUniformDistribution : public BaseDistribution<double>
    {
    public:
        using UniformDistributionType = std::uniform_real_distribution<double>;

        RealUniformDistribution()
        {
            m_LowerBoundParameter = 0.0;
            m_UpperBoundParameter = 1.0;
        }

        bool BelongsToSupport(const ValueType &x);
        double GetDensity(const ValueType &x);
        double GetLogDensity(const ValueType &x);
        double GetCumulative(const ValueType &x);
        void Fit(const SampleType &sample, const std::string &method);
        void Random(SampleType &sample, GeneratorType &generator);
        ValueType GetMean() { return (m_LowerBoundParameter + m_UpperBoundParameter) / 2.0; }
        double GetVariance() { return (m_UpperBoundParameter - m_LowerBoundParameter) * (m_UpperBoundParameter - m_LowerBoundParameter) / 12.0; }
        double GetDistance(Self *otherDistribution);

        void SetLowerBoundParameter(const double val) { m_LowerBoundParameter = val; }
        double GetLowerBoundParameter() { return m_LowerBoundParameter; }

        void SetUpperBoundParameter(const double val) { m_UpperBoundParameter = val; }
        double GetUpperBoundParameter() { return m_UpperBoundParameter; }

    private:
        double m_LowerBoundParameter;
        double m_UpperBoundParameter;
    };

} // end of namespace
