#pragma once

#include <vector>
#include <itkVariableLengthVector.h>

namespace anima
{

template <class ScalarType>
double
ComputeAsymmetricPDF(const std::vector<ScalarType> &w, const std::vector<ScalarType> &mu,
                     const ScalarType &kappa, const ScalarType &d, const ScalarType &nu,
                     const ScalarType &step, std::vector<ScalarType> &integrand);

template <class VectorType, class ScalarType>
double ComputeSymmetricPDF(const VectorType &w, const VectorType &mu,
                           const ScalarType &kappa, const ScalarType &d, const ScalarType &nu,
                           const ScalarType &step, std::vector<ScalarType> &integrand);

template <class T1, class VectorType>
double ComputeSymmetricCDF(const VectorType &direction, const T1 &kappa, const T1 &lambda,
                           const T1 &nu, const T1 &bvalue, const VectorType &gradient);

template <class T1, class T2, class VectorType>
double ComputeMixtureSymmetricCDF(const unsigned int &NbComponents,
                                  const T2 *Directions, const T1 *Kappa,
                                  const T1 *Lambda, const T1 *Nu, const T1 *W,
                                  const T1 &b, const VectorType &grad);

} // end of namespace anima

#include "animaDDIDistribution.hxx"
