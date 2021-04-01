#include <animaB1GMMDistributionIntegrand.h>

namespace anima
{

std::vector <double> B1GMMDistributionIntegrand::operator() (double const t)
{
    std::vector <double> epgVector = m_EPGSimulator.GetValue(m_T1Value, t, m_FlipAngle, 1.0);

    double gaussianExponent = (t - m_GaussianMean) * (t - m_GaussianMean) / (2.0 * m_GaussianVariance);

    for (unsigned int i = 0;i < epgVector.size();++i)
        epgVector[i] *= std::exp(- gaussianExponent) / (std::sqrt(2.0 * M_PI * m_GaussianVariance));

    return epgVector;
}

} // end namespace anima
