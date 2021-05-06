#include "animaB1GammaDistributionIntegrand.h"

#include <boost/math/distributions/gamma.hpp>
#include <boost/math/special_functions/polygamma.hpp>

namespace anima
{

std::vector <double> B1GammaDistributionIntegrand::operator() (double const t)
{
    std::vector <double> epgVector = m_EPGSimulator.GetValue(m_T1Value, t, m_FlipAngle, 1.0);

    double shape = m_GammaMean * m_GammaMean / m_GammaVariance;
    double scale = m_GammaVariance / m_GammaMean;

    double gammaValue = boost::math::gamma_p_derivative(shape, t / scale) / scale;

    for (unsigned int i = 0;i < epgVector.size();++i)
        epgVector[i] *= gammaValue;

    return epgVector;
}

} // end namespace anima
