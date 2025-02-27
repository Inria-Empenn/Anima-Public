#pragma once

#include <animaBaseDistribution.h>

#include <itkMatrix.h>
#include <itkVector.h>

namespace anima
{
    class ANIMASTATISTICALDISTRIBUTIONS_EXPORT WatsonDistribution : public BaseDistribution<itk::Vector<double, 3>>
    {
    public:
        using UniformDistributionType = std::uniform_real_distribution<double>;

        WatsonDistribution()
        {
            m_MeanAxis[0] = 0;
            m_MeanAxis[1] = 0;
            m_MeanAxis[2] = 1;
            m_ConcentrationParameter = 1.0;
            m_RValue = 0.0;
        }

        bool BelongsToSupport(const ValueType &x);
        double GetDensity(const ValueType &x);
        double GetLogDensity(const ValueType &x);
        double GetCumulative(const ValueType &x);
        void Fit(const SampleType &sample, const std::string &method);
        void Random(SampleType &sample, GeneratorType &generator);
        ValueType GetMean();
        double GetVariance() { return 1.0 - m_RValue; }
        double GetDistance(Self *otherDistribution);

        void SetMeanAxis(const ValueType &x);
        ValueType GetMeanAxis() { return m_MeanAxis; }

        void SetConcentrationParameter(const double &x);
        double GetConcentrationParameter() { return m_ConcentrationParameter; }

        vnl_matrix<double> GetCovarianceMatrix();
        void GetStandardWatsonSHCoefficients(
            std::vector<double> &coefficients,
            std::vector<double> &derivatives);

    private:
        double ComputeConcentrationMLE(const double rValue, const double aValue, const double cValue, double &logLik);
        ValueType m_MeanAxis;
        double m_ConcentrationParameter;
        double m_RValue;
        itk::Matrix<double, 3, 3> m_NorthToMeanAxisRotationMatrix;
        const unsigned int m_AmbientDimension = 3;
    };

} // end of namespace anima
