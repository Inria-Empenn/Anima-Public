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

    for (int k = 0;k <= (int)m_LOrder;k += 2)
        for (int m = -k;m <= k;++m)
        {
            SphericalHarmonic tmpSH(k,m);
            m_SphericalHarmonics.push_back(tmpSH);
        }
}

double ODFSphericalHarmonicBasis::getNthSHValueAtPosition(int k, int m, double theta, double phi)
{
    std::complex <double> tmpVal = m_SphericalHarmonics[k*(k+1)/2 + m].Value(theta,phi);

    double resVal = 0;

    if (m > 0)
        resVal = sqrt(2.0)*imag(tmpVal);
    else if (m < 0)
        resVal = sqrt(2.0)*real(tmpVal);
    else
        resVal = real(tmpVal);

    return resVal;
}

} // end namespace anima
