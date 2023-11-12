#pragma once

#include <animaBaseDistribution.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_trace.h>

namespace anima
{
    class ANIMASTATISTICALDISTRIBUTIONS_EXPORT MultivariateNormalDistribution : public BaseDistribution<std::vector<double>>
    {
    public:
        using NormalDistributionType = std::normal_distribution<double>;
        using MatrixType = vnl_matrix<double>;

        MultivariateNormalDistribution()
        {
            m_MeanParameter.clear();
            m_CovarianceMatrixParameter.clear();
            m_SqrtCovarianceMatrixParameter.clear();
            m_PrecisionMatrixParameter.clear();
            m_CovarianceMatrixDeterminant = 1.0;
        }

        bool BelongsToSupport(const ValueType &x) { return true; }
        double GetDensity(const ValueType &x);
        double GetLogDensity(const ValueType &x);
        double GetCumulative(const ValueType &x);
        void Fit(const SampleType &sample, const std::string &method);
        void Random(SampleType &sample, GeneratorType &generator);
        ValueType GetMean() { return m_MeanParameter; }
        double GetVariance() { return vnl_trace(this->GetCovarianceMatrixParameter()); }
        double GetDistance(Self *otherDistribution);

        void SetMeanParameter(const ValueType &val) { m_MeanParameter = val; }
        ValueType GetMeanParameter() { return m_MeanParameter; }

        void SetCovarianceMatrixParameter(const MatrixType &val);
        MatrixType GetCovarianceMatrixParameter() { return m_CovarianceMatrixParameter; }

        MatrixType GetPrecisionMatrixParameter() { return m_PrecisionMatrixParameter; }

    private:
        ValueType m_MeanParameter;
        MatrixType m_CovarianceMatrixParameter, m_SqrtCovarianceMatrixParameter, m_PrecisionMatrixParameter;
        double m_CovarianceMatrixDeterminant;
    };

} // end of namespace
