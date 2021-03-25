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

    void SetL(int &l) { if (m_L != l) {m_L = l; m_NeedUpdate = true;} }
    void SetM(int &m) { if (m_M != m) {m_M = m; m_NeedUpdate = true;} }

    std::complex <double> Value(const double &theta, const double &phi);

    std::complex <double> getThetaFirstDerivative(const double& theta, const double& phi);
    std::complex <double> getPhiFirstDerivative(const double& theta, const double& phi);

    std::complex <double> getThetaSecondDerivative(const double& theta, const double& phi);
    std::complex <double> getPhiSecondDerivative(const double& theta, const double& phi);

    std::complex <double> getThetaPhiDerivative(const double& theta, const double& phi);

private:
    void UpdateSQRTFactor();
    int m_L;
    int m_M;

    double m_SQRTFactor;
    bool m_NeedUpdate;
};

} // end of namespace anima


