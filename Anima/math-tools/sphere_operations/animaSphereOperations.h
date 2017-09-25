#pragma once

#include <vector>

namespace anima
{

/**
 * Computes N samples evenly spread on the sphere using the Fibonacci spiral algorithm.
 *
 * @param spherePoints output vector containing all sampled point
 * @param numSamples number of output samples required
 */
template <class ScalarType> void GetSphereEvenSampling(std::vector < std::vector <ScalarType> > &spherePoints,
                                                       unsigned int numSamples);

} // end namespace anima

#include "animaSphereOperations.hxx"
