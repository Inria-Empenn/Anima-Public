#pragma once
#include "animaODFSphericalHarmonicBasis.h"

namespace anima
{

template <class T>
double
ODFSphericalHarmonicBasis::
getValueAtPosition(const T &coefficients, double theta, double phi)
{
    double resVal = 0;
    std::complex <double> tmpVal;

    for (int k = 0;k <= (int)m_LOrder;k += 2)
        for (int m = -k;m <= k;++m)
        {
            tmpVal = m_SphericalHarmonics[k*(k+1)/2 + m].Value(theta,phi);

            if (m > 0)
                resVal += coefficients[k*(k+1)/2 + m]*sqrt(2.0)*imag(tmpVal);
            else if (m < 0)
                resVal += coefficients[k*(k+1)/2 + m]*sqrt(2.0)*real(tmpVal);
            else
                resVal += coefficients[k*(k+1)/2 + m]*real(tmpVal);
        }

    return resVal;
}

template <class T>
itk::VariableLengthVector <T>
ODFSphericalHarmonicBasis::
GetSampleValues(itk::VariableLengthVector <T> &data,
                std::vector < std::vector <double> > &m_SampleDirections)
{
    itk::VariableLengthVector <T> resVal(m_SampleDirections.size());

    for (unsigned int i = 0;i < m_SampleDirections.size();++i)
        resVal[i] = this->getValueAtPosition(data,m_SampleDirections[i][0],m_SampleDirections[i][1]);

    return resVal;
}

template <class T>
double
ODFSphericalHarmonicBasis::
getThetaFirstDerivativeValueAtPosition(const T &coefficients, double theta, double phi)
{
    double resVal = 0;
    std::complex <double> tmpVal;

    for (int k = 0;k <= (int)m_LOrder;k += 2)
        for (int m = -k;m <= k;++m)
        {
            tmpVal = m_SphericalHarmonics[k*(k+1)/2 + m].getThetaFirstDerivative(theta,phi);

            if (m > 0)
                resVal += coefficients[k*(k+1)/2 + m]*sqrt(2.0)*imag(tmpVal);
            else if (m < 0)
                resVal += coefficients[k*(k+1)/2 + m]*sqrt(2.0)*real(tmpVal);
            else
                resVal += coefficients[k*(k+1)/2 + m]*real(tmpVal);
        }

    return resVal;
}

template <class T>
double
ODFSphericalHarmonicBasis::
getPhiFirstDerivativeValueAtPosition(const T &coefficients, double theta, double phi)
{
    double resVal = 0;
    std::complex <double> tmpVal;

    for (int k = 0;k <= (int)m_LOrder;k += 2)
        for (int m = -k;m <= k;++m)
        {
            tmpVal = m_SphericalHarmonics[k*(k+1)/2 + m].getPhiFirstDerivative(theta,phi);

            if (m > 0)
                resVal += coefficients[k*(k+1)/2 + m]*sqrt(2.0)*imag(tmpVal);
            else if (m < 0)
                resVal += coefficients[k*(k+1)/2 + m]*sqrt(2.0)*real(tmpVal);
            else
                resVal += coefficients[k*(k+1)/2 + m]*real(tmpVal);
        }

    return resVal;
}

template <class T>
double
ODFSphericalHarmonicBasis::
getThetaSecondDerivativeValueAtPosition(const T &coefficients, double theta, double phi)
{
    double resVal = 0;
    std::complex <double> tmpVal;

    for (int k = 0;k <= (int)m_LOrder;k += 2)
        for (int m = -k;m <= k;++m)
        {
            tmpVal = m_SphericalHarmonics[k*(k+1)/2 + m].getThetaSecondDerivative(theta,phi);

            if (m > 0)
                resVal += coefficients[k*(k+1)/2 + m]*sqrt(2.0)*imag(tmpVal);
            else if (m < 0)
                resVal += coefficients[k*(k+1)/2 + m]*sqrt(2.0)*real(tmpVal);
            else
                resVal += coefficients[k*(k+1)/2 + m]*real(tmpVal);
        }

    return resVal;
}

template <class T>
double
ODFSphericalHarmonicBasis::
getThetaPhiDerivativeValueAtPosition(const T &coefficients, double theta, double phi)
{
    double resVal = 0;
    std::complex <double> tmpVal;

    for (int k = 0;k <= (int)m_LOrder;k += 2)
        for (int m = -k;m <= k;++m)
        {
            tmpVal = m_SphericalHarmonics[k*(k+1)/2 + m].getThetaPhiDerivative(theta,phi);

            if (m > 0)
                resVal += coefficients[k*(k+1)/2 + m]*sqrt(2.0)*imag(tmpVal);
            else if (m < 0)
                resVal += coefficients[k*(k+1)/2 + m]*sqrt(2.0)*real(tmpVal);
            else
                resVal += coefficients[k*(k+1)/2 + m]*real(tmpVal);
        }

    return resVal;
}

template <class T>
double
ODFSphericalHarmonicBasis::
getPhiSecondDerivativeValueAtPosition(const T &coefficients, double theta, double phi)
{
    double resVal = 0;
    std::complex <double> tmpVal;

    for (int k = 0;k <= (int)m_LOrder;k += 2)
        for (int m = -k;m <= k;++m)
        {
            tmpVal = m_SphericalHarmonics[k*(k+1)/2 + m].getPhiSecondDerivative(theta,phi);

            if (m > 0)
                resVal += coefficients[k*(k+1)/2 + m]*sqrt(2.0)*imag(tmpVal);
            else if (m < 0)
                resVal += coefficients[k*(k+1)/2 + m]*sqrt(2.0)*real(tmpVal);
            else
                resVal += coefficients[k*(k+1)/2 + m]*real(tmpVal);
        }

    return resVal;
}

template <class T>
double
ODFSphericalHarmonicBasis::
getCurvatureAtPosition(const T &coefficients, double theta, double phi)
{
    // Taken from Bloy and Verma, simplified to the maximum (this supposes we are actually at an extremum of the odf)
    double odfValue = this->getValueAtPosition(coefficients,theta,phi);
    double sqSinTheta = sin(theta) * sin(theta);
    if (sqSinTheta <= 1.0e-16)
        sqSinTheta = 1.0e-16;

    double denom = 2.0 * odfValue * odfValue * sqSinTheta;

    double num = 2.0 * odfValue * sqSinTheta;
    num -= sqSinTheta * this->getThetaSecondDerivativeValueAtPosition(coefficients,theta,phi) +
    this->getPhiSecondDerivativeValueAtPosition(coefficients,theta,phi);

    return num / denom;
}

} // end of namespace anima
