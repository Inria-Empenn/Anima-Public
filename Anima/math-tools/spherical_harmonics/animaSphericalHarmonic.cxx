#include <cmath>

#include "animaSphericalHarmonic.h"

#include <boost/math/special_functions/legendre.hpp>
#include <animaLegendreDerivatives.h>

namespace anima
{

SphericalHarmonic::SphericalHarmonic()
{
    m_L = 0;
    m_M = 0;

    m_NeedUpdate = true;
}

SphericalHarmonic::SphericalHarmonic(int &l, int &m)
{
    m_L = l;
    m_M = m;

    m_NeedUpdate = true;
}

void SphericalHarmonic::UpdateSQRTFactor()
{
    if (!m_NeedUpdate)
        return;

    int absm = std::abs(m_M);
    double factorialRatio = 1;
    for (int i = m_L - absm + 1;i <= m_L + absm;++i)
        factorialRatio *= i;

    m_SQRTFactor = std::sqrt((2*m_L + 1) / (4 * M_PI * factorialRatio));
    m_NeedUpdate = false;
}

std::complex <double> SphericalHarmonic::Value(const double &theta, const double &phi)
{
    int absm = std::abs(m_M);
    std::complex <double> resVal(0, absm * phi);
    resVal = std::exp(resVal);

    this->UpdateSQRTFactor();

    resVal *= m_SQRTFactor * boost::math::legendre_p(m_L,absm,std::cos(theta));
    if (m_M < 0)
    {
        resVal = std::conj(resVal);

        if (absm % 2 != 0)
            resVal *= -1;
    }

    return resVal;
}

std::complex <double> SphericalHarmonic::getThetaFirstDerivative(const double& theta, const double& phi)
{
    if ((m_M == 0)||(m_L < 1))
        return std::complex<double> (0.0,0.0);

    int absm = std::abs(m_M);
    std::complex<double> retval(0.0,(double)(absm*phi));
    retval = std::exp(retval);

    this->UpdateSQRTFactor();

    double factor = m_SQRTFactor;
    if ((absm%2 != 0)&&(m_M < 0))
        factor *= -1;

    factor *= - sin(theta) * anima::legendre_first_derivative(m_L,absm,cos(theta));

    return factor * retval;
}

std::complex <double> SphericalHarmonic::getPhiFirstDerivative(const double& theta, const double& phi)
{
    int absm = std::abs(m_M);
    std::complex<double> retval(0.0,(double)(absm*phi));
    retval = std::exp(retval);
    retval *= std::complex<double> (0.0,absm);

    this->UpdateSQRTFactor();

    double factor = m_SQRTFactor * boost::math::legendre_p (m_L,absm,cos(theta));

    if ((absm%2 != 0)&&(m_M < 0))
        factor *= -1;

    return factor * retval;
}

std::complex <double> SphericalHarmonic::getThetaSecondDerivative(const double& theta, const double& phi)
{
    if ((m_M == 0)||(m_L < 2))
        return std::complex<double> (0.0,0.0);

    int absm = std::abs(m_M);
    std::complex<double> retval(0.0,(double)(absm*phi));
    retval = std::exp(retval);

    this->UpdateSQRTFactor();

    double factor = m_SQRTFactor;
    if ((absm%2 != 0)&&(m_M < 0))
        factor *= -1;

    factor *= sin(theta) * sin(theta) * anima::legendre_second_derivative(m_L,absm,cos(theta)) -
    cos(theta) * anima::legendre_first_derivative(m_L,absm,cos(theta));

    return factor * retval;
}

std::complex <double> SphericalHarmonic::getPhiSecondDerivative(const double& theta, const double& phi)
{
    double factor = - m_M * m_M;

    return factor * this->Value(theta,phi);
}

std::complex <double> SphericalHarmonic::getThetaPhiDerivative(const double& theta, const double& phi)
{
    if ((m_M == 0)||(m_L < 1))
        return std::complex<double> (0.0,0.0);

    int absm = std::abs(m_M);
    std::complex<double> retval(0.0,(double)(absm*phi));
    retval = std::exp(retval);
    retval *= std::complex<double> (0.0,absm);

    this->UpdateSQRTFactor();

    double factor = m_SQRTFactor;
    if ((absm%2 != 0)&&(m_M < 0))
        factor *= -1;

    factor *= sin(theta) * anima::legendre_first_derivative(m_L,absm,cos(theta));

    return factor * retval;
}

} // end of namespace anima
