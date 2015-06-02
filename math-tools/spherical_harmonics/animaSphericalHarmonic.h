#pragma once

#include <complex>
#include "AnimaSHToolsExport.h"

namespace anima
{

class ANIMASHTOOLS_EXPORT SphericalHarmonic
{
public:
    SphericalHarmonic();
    SphericalHarmonic(int &l, int &m);

    void SetL(int &l) {m_L = l;}
    void SetM(int &m) {m_M = m;}

    std::complex <double> Value(const double &theta, const double &phi);

    std::complex <double> getThetaFirstDerivative(const double& theta, const double& phi);
    std::complex <double> getPhiFirstDerivative(const double& theta, const double& phi);

    std::complex <double> getThetaSecondDerivative(const double& theta, const double& phi);
    std::complex <double> getPhiSecondDerivative(const double& theta, const double& phi);

    std::complex <double> getThetaPhiDerivative(const double& theta, const double& phi);

private:
    int m_L;
    int m_M;
};

} // end of namespace anima


