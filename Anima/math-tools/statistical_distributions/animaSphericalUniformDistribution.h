#pragma once

#include <animaBaseDistribution.h>

#include <itkVector.h>

namespace anima
{
    class ANIMASTATISTICALDISTRIBUTIONS_EXPORT SphericalUniformDistribution : public BaseDistribution<itk::Vector<double, 3>>
    {
    public:
        using UniformDistributionType = std::uniform_real_distribution<double>;

        SphericalUniformDistribution() {}

        bool BelongsToSupport(const ValueType &x);
        double GetDensity(const ValueType &x);
        double GetLogDensity(const ValueType &x);
        double GetCumulative(const ValueType &x);
        void Fit(const SampleType &sample, const std::string &method) { return; }
        void Random(SampleType &sample, GeneratorType &generator);
        ValueType GetMean();
        double GetVariance() { return 0.0; }
        double GetDistance(Self *otherDistribution);
    };

} // end of namespace
