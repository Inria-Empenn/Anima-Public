#include "animaB1GammaDerivativeDistributionIntegrand.h"

#include <boost/math/special_functions/polygamma.hpp>

namespace anima
{

double B1GammaDerivativeDistributionIntegrand::operator() (double const t)
{
    if ((m_EPGVectors.find(t) == m_EPGVectors.end())||(m_B1DerivativeFlag && (m_EPGVectors.find(t) == m_EPGVectors.end())))
    {
        m_EPGVectors.insert(std::make_pair(t,m_EPGSimulator.GetValue(m_T1Value, t, m_FlipAngle, 1.0)));
        if (m_B1DerivativeFlag)
            m_DerivativeEPGVectors.insert(std::make_pair(t,m_EPGSimulator.GetFADerivative()));
    }

    if (m_B1DerivativeFlag)
    {
        // Derivative against flip angle parameter
        double shape = m_GammaMean * m_GammaMean / m_GammaVariance;
        double scale = m_GammaVariance / m_GammaMean;

        double gammaValue = boost::math::gamma_p_derivative(shape, t / scale) / scale;

        return m_DerivativeEPGVectors[t][m_EchoNumber] * gammaValue;
    }
    else
    {
        // Derivative against mean parameter of gamma distribution
        double shape = m_GammaMean * m_GammaMean / m_GammaVariance;
        double scale = m_GammaVariance / m_GammaMean;

        double internalTerm = (2.0 * std::log(t / scale) - 2.0 * boost::math::digamma(shape) + 1.0) / scale - t / m_GammaVariance;
        double derivativeGammaValue = internalTerm * boost::math::gamma_p_derivative(shape, t / scale) / scale;

        return m_EPGVectors[t][m_EchoNumber] * derivativeGammaValue;
    }
}

} // end namespace anima
