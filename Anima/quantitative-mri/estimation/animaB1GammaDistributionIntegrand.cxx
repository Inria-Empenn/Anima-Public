#include "animaB1GammaDistributionIntegrand.h"

#include <boost/math/distributions/gamma.hpp>
#include <boost/math/special_functions/polygamma.hpp>

namespace anima
{

double B1GammaDistributionIntegrand::operator() (double const t)
{
    if (m_EPGVectors.find(t) == m_EPGVectors.end())
        m_EPGVectors.insert(std::make_pair(t,m_EPGSimulator.GetValue(m_T1Value, t, m_FlipAngle, 1.0)));

    double shape = m_GammaMean * m_GammaMean / m_GammaVariance;
    double scale = m_GammaVariance / m_GammaMean;

    double gammaValue = boost::math::gamma_p_derivative(shape, t / scale) / scale;

    return m_EPGVectors[t][m_EchoNumber] * gammaValue;
}

} // end namespace anima
