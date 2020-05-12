#include <animaB1GMMDistributionIntegrand.h>

namespace anima
{

double B1GMMDistributionIntegrand::operator() (double const t)
{
    if (m_EPGVectors.find(t) == m_EPGVectors.end())
        m_EPGVectors.insert(std::make_pair(t,m_EPGSimulator.GetValue(m_T1Value, t, m_FlipAngle, 1.0)));

    double gaussianExponent = (t - m_GaussianMean) * (t - m_GaussianMean) / (2.0 * m_GaussianVariance);

    return m_EPGVectors[t][m_EchoNumber] * std::exp(- gaussianExponent) / (std::sqrt(2.0 * M_PI * m_GaussianVariance));
}

} // end namespace anima
