#pragma once

#include <vector>

namespace anima
{

/**
 * Computes N samples almost evenly spread on the half-sphere.
 * C.G. Koay. A simple scheme for generating nearly uniform distribution of antipodally symmetric pointson the unit sphere
 * Journal of computer science, 2, 2011, pp. 377-381
 *
 * @param spherePoints output vector containing all sampled point
 * @param numSamples number of output samples required
 */
template <class ScalarType> void GetSphereEvenSampling(std::vector < std::vector <ScalarType> > &spherePoints,
                                                       unsigned int numSamples);

} // end namespace anima

#include "animaSphereOperations.hxx"
