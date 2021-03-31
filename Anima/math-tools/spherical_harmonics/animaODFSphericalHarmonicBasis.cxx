#include "animaODFSphericalHarmonicBasis.h"

namespace anima
{

ODFSphericalHarmonicBasis::ODFSphericalHarmonicBasis(unsigned int L)
{
    m_LOrder = 0;
    this->SetOrder(L);
}

void ODFSphericalHarmonicBasis::SetOrder(unsigned int L)
{
    if (m_LOrder == L)
        return;

    m_LOrder = L;
    m_SphericalHarmonics.clear();

    for (unsigned int k = 0;k <= m_LOrder;k += 2)
        for (unsigned int m = 0;m <= k;++m)
        {
            SphericalHarmonic tmpSH(k,m);
            m_SphericalHarmonics.push_back(tmpSH);
        }
}

double ODFSphericalHarmonicBasis::getNthSHValueAtPosition(int k, int m, double theta, double phi)
{
    int absm = std::abs(m);
    std::complex <double> tmpVal = m_SphericalHarmonics[k * k / 4 + absm].Value(theta,phi);

    double resVal = 0;

    if (m > 0)
        resVal = sqrt(2.0)*imag(tmpVal);
    else if (m < 0)
    {
        resVal = sqrt(2.0)*real(tmpVal);
        if (m % 2 != 0)
            resVal *= -1;
    }
    else
        resVal = real(tmpVal);

    return resVal;
}

} // end namespace anima
