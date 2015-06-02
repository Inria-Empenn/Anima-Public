#pragma once

#include <vector>

namespace anima
{

template <class T> T Cubic (T x0, T x1, T x2, T x3, T y0, T y1, T y2, T y3, T x);
template <class T> void InverseCubicInterpolator(std::vector <T> &inputVect, std::vector <T> &outputVect, T step);
template <class T> void CubicInterpolator(std::vector <T> &transfVect, std::vector <T> &scale, std::vector<T> &outputVect, T LengthLine);

} // end namespace anima

#include "animaCubicInterpolation.hxx"
