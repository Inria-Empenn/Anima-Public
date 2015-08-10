#pragma once

#include <boost/random/mersenne_twister.hpp>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector_fixed.h>

#include <itkPoint.h>
#include <itkVector.h>

namespace anima
{

template <class T>
double SampleFromUniformDistribution(const T &a, const T &b, boost::mt19937 &generator);

template <class T>
unsigned int SampleFromBernoulliDistribution(const T &p, boost::mt19937 &generator);

template <class T>
double SampleFromGaussianDistribution(const T &mean, const T &std, boost::mt19937 &generator);

template <class VectorType, class ScalarType>
void SampleFromMultivariateGaussianDistribution(const VectorType &mean, const vnl_matrix <ScalarType> &mat, VectorType &resVec,
                                                boost::mt19937 &generator, bool isMatCovariance = true);

// From Ulrich 1984
template <class VectorType, class ScalarType>
void SampleFromVMFDistribution(const ScalarType &kappa, const VectorType &meanDirection, VectorType &resVec, boost::mt19937 &generator);

// From Wenzel 2012
template <class VectorType, class ScalarType>
void SampleFromVMFDistributionNumericallyStable(const ScalarType &kappa, const VectorType &meanDirection, VectorType &resVec, boost::mt19937 &generator);

template <class ScalarType, class VectorType>
void
SampleFromWatsonDistribution(const ScalarType &kappa, const VectorType &meanDirection, VectorType &resVec, unsigned int DataDimension, boost::mt19937 &generator);

template <class ScalarType, unsigned int DataDimension>
void
SampleFromWatsonDistribution(const ScalarType &kappa, const vnl_vector_fixed < ScalarType, DataDimension > &meanDirection, vnl_vector_fixed < ScalarType, DataDimension > &resVec, boost::mt19937 &generator);

template <class ScalarType, unsigned int DataDimension>
void
SampleFromWatsonDistribution(const ScalarType &kappa, const itk::Point < ScalarType, DataDimension > &meanDirection, itk::Point < ScalarType, DataDimension > &resVec, boost::mt19937 &generator);

template <class ScalarType, unsigned int DataDimension>
void
SampleFromWatsonDistribution(const ScalarType &kappa, const itk::Vector < ScalarType, DataDimension > &meanDirection, itk::Vector < ScalarType, DataDimension > &resVec, boost::mt19937 &generator);

} // end of namespace anima

#include "animaDistributionSampling.hxx"
