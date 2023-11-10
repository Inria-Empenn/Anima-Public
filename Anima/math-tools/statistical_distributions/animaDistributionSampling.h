#pragma once

#include <random>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector_fixed.h>

#include <itkPoint.h>
#include <itkVector.h>

namespace anima
{

template <class VectorType>
void SampleFromUniformDistributionOn2Sphere(std::mt19937 &generator, VectorType &resVec);

template <class T>
unsigned int SampleFromBernoulliDistribution(const T &p, std::mt19937 &generator);

template <class T>
double SampleFromGaussianDistribution(const T &mean, const T &std, std::mt19937 &generator);

template <class VectorType, class ScalarType>
void SampleFromMultivariateGaussianDistribution(const VectorType &mean, const vnl_matrix <ScalarType> &mat, VectorType &resVec,
                                                std::mt19937 &generator, bool isMatCovariance = true);

} // end of namespace anima

#include "animaDistributionSampling.hxx"
