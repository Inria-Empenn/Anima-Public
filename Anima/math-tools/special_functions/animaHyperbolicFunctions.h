#pragma once

namespace anima
{

template <class T> double ShOverId(const T &x);
template <class ScalarType> std::complex<double> ComplexShRatio(const ScalarType &k, const ScalarType &alpha, const ScalarType &beta);
template <class T1, class T2, class T3> double ShRatio(const T1 &k, const T2 &alpha, const T3 &beta);
template <class T> double xi(const T &k);
template <class T> double jtwo(const T &k);
template <class T> double jfour(const T &k);

} // end namespace hyperbolic_functions

#include "animaHyperbolicFunctions.hxx"


