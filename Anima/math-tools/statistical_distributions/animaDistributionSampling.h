#pragma once

#include <random>
#include <vnl/vnl_matrix.h>

namespace anima
{

    template <class VectorType>
    void SampleFromUniformDistributionOn2Sphere(std::mt19937 &generator, VectorType &resVec);

    template <class VectorType, class ScalarType>
    void SampleFromMultivariateGaussianDistribution(const VectorType &mean, const vnl_matrix<ScalarType> &mat, VectorType &resVec,
                                                    std::mt19937 &generator, bool isMatCovariance = true);

} // end of namespace anima

#include "animaDistributionSampling.hxx"
