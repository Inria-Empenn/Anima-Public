#include <animaGaussLaguerreQuadrature.h>
#include <cmath>
#include <iostream>

#include <boost/math/distributions/gamma.hpp>
#include <boost/math/quadrature/gauss.hpp>

struct IntegrandType
{
    double operator()(const double t)
    {
        double shape = m_GammaMean * m_GammaMean / m_GammaVariance;
        double scale = m_GammaVariance / m_GammaMean;

        double gammaValue = boost::math::gamma_p_derivative(shape, t / scale) / scale;

        return gammaValue;
    }

    double m_GammaMean;
    double m_GammaVariance;
};

int main(int argc, char **argv)
{
    double gammaMean = 15.0;
    double gammaVariance = 50.0;

    anima::GaussLaguerreQuadrature glQuad;

    IntegrandType integrand;
    integrand.m_GammaMean = gammaMean;
    integrand.m_GammaVariance = gammaVariance;

    for (unsigned int i = 1;i < 11;++i)
    {
        double minValue = std::max(0.0,gammaMean - i * std::sqrt(gammaVariance));
        double theoreticalMinValue = gammaMean - i * std::sqrt(gammaVariance);
        double maxValue = gammaMean + i * std::sqrt(gammaVariance);
        if (theoreticalMinValue < 0.0)
            maxValue -= theoreticalMinValue;
        std::cout << "Interest zone " << i << " [" << minValue << "," << maxValue << "]" << std::endl;
        glQuad.SetInterestZone(minValue, maxValue);

        double integralValue = glQuad.GetIntegralValue(integrand);
        std::cout << "Gauss Laguerre value: " << integralValue << std::endl;

        integralValue = boost::math::quadrature::gauss <double, 15>::integrate(integrand, minValue, maxValue);
        std::cout << "Gauss Legendre value: " << integralValue << std::endl;
    }

    return EXIT_SUCCESS;
}
