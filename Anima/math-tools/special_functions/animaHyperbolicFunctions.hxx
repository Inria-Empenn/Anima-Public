#pragma once

#include "animaHyperbolicFunctions.h"

namespace anima
{

template <class T> double ShOverId(const T &x)
{
    double resVal = x;

    if (resVal < 1e-4)
        resVal = 1.0 + resVal * resVal / 6.0;
    else
        resVal = std::sinh(resVal) / resVal;

    return resVal;
}

template <class ScalarType> std::complex<double> ComplexShRatio(const ScalarType &k, const ScalarType &alpha, const ScalarType &beta)
{
    double ss = alpha * alpha + beta * beta;
    
    double constant = k * std::exp(alpha - k) / ss / (1 - std::exp(- 2.0 * k));
    
    double realVal = constant * (alpha * (1.0 - std::exp(- 2.0 * alpha)) * std::cos(beta) + beta * (1.0 + std::exp(- 2.0 * alpha)) * std::sin(beta));
    
    double imagVal = constant * (alpha * (1.0 + std::exp(- 2.0 * alpha)) * std::sin(beta) - beta * (1.0 - std::exp(- 2.0 * alpha)) * std::cos(beta));
    
    return std::complex<double>(realVal, imagVal);
}
    
template <class T1, class T2, class T3> double ShRatio(const T1 &k, const T2 &alpha, const T3 &beta)
{
    double ss = alpha * alpha + beta * beta;

    double resVal = k * exp(alpha - k) * (alpha * (1.0 - exp(- 2.0 * alpha)) * cos(beta) + beta * (1.0 + exp(- 2.0 * alpha)) * sin(beta)) / ss / (1 - exp(- 2.0 * k));

    return resVal;
}

template <class T> double  xi(const T &k)
{
    double resVal = k * k;

    if (k < 1e-4)
        resVal = 1.0 / 3.0 - resVal / 45;
    else
        resVal = 1.0 / (k * tanh(k)) - 1.0 / resVal;

    return resVal;
}

template <class T> double jtwo(const T &k)
{
    double jtwoVal = k * k;

    if (k < 1e-4)
        jtwoVal = 1.0 / 3.0 + 2.0 * k * k / 45.0;
    else
        jtwoVal = 1.0 - 2.0 * (1.0 / (k * tanh(k)) - 1.0 / jtwoVal);

    return jtwoVal;
}

template <class T> double jfour(const T &k)
{
    double jfourVal = k * k;

    if (k < 1e-3)
        jfourVal = 1.0 / 5.0 + 4.0 * k * k / 45.0;
    else
        jfourVal = 1.0 - 4.0 * (1.0 / (k * tanh(k)) - 3.0 * (1.0 - 2.0 * (1.0 / (k * tanh(k)) - 1.0 / jfourVal)) / jfourVal);

    return jfourVal;
}

} // end of namespace anima


