#pragma once

#include <random>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector_fixed.h>

#include <itkPoint.h>
#include <itkVector.h>

namespace anima
{

template <class T>
double SampleFromUniformDistribution(const T &a, const T &b, std::mt19937 &generator);

template <class VectorType>
void SampleFromUniformDistributionOn2Sphere(std::mt19937 &generator, VectorType &resVec);

template <class T>
unsigned int SampleFromBernoulliDistribution(const T &p, std::mt19937 &generator);

template <class T>
double SampleFromGaussianDistribution(const T &mean, const T &std, std::mt19937 &generator);

template <class VectorType, class ScalarType>
void SampleFromMultivariateGaussianDistribution(const VectorType &mean, const vnl_matrix <ScalarType> &mat, VectorType &resVec,
                                                std::mt19937 &generator, bool isMatCovariance = true);

// From Ulrich 1984
template <class VectorType, class ScalarType>
void SampleFromVMFDistribution(const ScalarType &kappa, const VectorType &meanDirection, VectorType &resVec, std::mt19937 &generator);

// From Wenzel 2012
template <class VectorType, class ScalarType>
void SampleFromVMFDistributionNumericallyStable(const ScalarType &kappa, const VectorType &meanDirection, VectorType &resVec, std::mt19937 &generator);

} // end of namespace anima

#include "animaDistributionSampling.hxx"
