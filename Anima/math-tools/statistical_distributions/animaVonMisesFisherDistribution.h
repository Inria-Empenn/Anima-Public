#pragma once

#include <animaBaseDistribution.h>

#include <itkVector.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_trace.h>

#include <boost/math/distributions/beta.hpp>

namespace anima
{
    class ANIMASTATISTICALDISTRIBUTIONS_EXPORT VonMisesFisherDistribution : public BaseDistribution<itk::Vector<double, 3>>
    {
    public:
        using UniformDistributionType = std::uniform_real_distribution<double>;
        using BetaDistributionType = boost::math::beta_distribution<double>;

        VonMisesFisherDistribution()
        {
            m_MeanDirection[0] = 0;
            m_MeanDirection[1] = 0;
            m_MeanDirection[2] = 1;
            m_ConcentrationParameter = 1.0;
            m_BesselRatio = 1.0;
        }

        bool BelongsToSupport(const ValueType &x);
        double GetDensity(const ValueType &x);
        double GetLogDensity(const ValueType &x);
        double GetCumulative(const ValueType &x);
        void Fit(const SampleType &sample, const std::string &method);
        void Random(SampleType &sample, GeneratorType &generator);
        ValueType GetMean() { return m_BesselRatio * m_MeanDirection; }
        double GetVariance() { return vnl_trace(this->GetCovarianceMatrix()); }
        double GetDistance(Self *otherDistribution);

        void SetMeanDirection(const ValueType &x);
        ValueType GetMeanDirection() { return m_MeanDirection; }

        void SetConcentrationParameter(const double &x);
        double GetConcentrationParameter() { return m_ConcentrationParameter; }

        double GetBesselRatio() { return m_BesselRatio; }
        vnl_matrix<double> GetCovarianceMatrix();

    private:
        void SampleFromVMFDistribution(ValueType &resVec, GeneratorType &generator);
        void SampleFromVMFDistributionNumericallyStable(ValueType &resVec, GeneratorType &generator);
        ValueType m_MeanDirection;
        double m_ConcentrationParameter;
        double m_BesselRatio;
        const unsigned int m_AmbientDimension = 3;
    };

} // end of namespace anima
