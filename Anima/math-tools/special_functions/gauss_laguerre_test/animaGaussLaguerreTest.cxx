#include <animaGaussLaguerreQuadrature.h>
#include <cmath>
#include <iostream>

#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/special_functions/erf.hpp>

struct IntegrandType
{
    double operator()(const double t)
    {
        double gaussianExponent = (t - m_GaussianMean) * (t - m_GaussianMean) / (2.0 * m_GaussianVariance);

        return std::exp(- gaussianExponent) / std::sqrt(2.0 * M_PI * m_GaussianVariance);
    }

    double m_GaussianMean;
    double m_GaussianVariance;
};

int main(int argc, char **argv)
{
    double gaussianMean = 20.0;
    double gaussianVariance = 25.0;

    anima::GaussLaguerreQuadrature glQuad;

    IntegrandType integrand;
    integrand.m_GaussianMean = gaussianMean;
    integrand.m_GaussianVariance = gaussianVariance;

    for (unsigned int i = 1;i < 11;++i)
    {
        double minValue = gaussianMean - i * std::sqrt(gaussianVariance);
        double maxValue = gaussianMean + i * std::sqrt(gaussianVariance);

        std::cout << "Interest zone " << i << " [" << minValue << "," << maxValue << "]" << std::endl;
        glQuad.SetInterestZone(minValue, maxValue);

        double integralValue = glQuad.GetIntegralValue(integrand);
        std::cout << "Gauss Laguerre value: " << integralValue << " " << boost::math::erf(gaussianMean / std::sqrt(2.0 * gaussianVariance))
                  << ", true value: " << 0.5 * (1.0 + boost::math::erf(gaussianMean / std::sqrt(2.0 * gaussianVariance))) << std::endl;

        integralValue = boost::math::quadrature::gauss <double, 15>::integrate(integrand, std::max(0.0,minValue), maxValue);
        std::cout << "Gauss Legendre value: " << integralValue << std::endl;
    }

    return EXIT_SUCCESS;
}
