#pragma once

#include <vector>

namespace anima
{

/**
 * Computes the logarithm map on the hyper-sphere.
 *
 * @param p point to compute the logarithm from
 * @param x reference point on the sphere (unfolding happens from that point)
 * @param logValue Logarithm value (result)
 */
template <class ScalarType> void sphere_log_map(const std::vector <ScalarType> &p, const std::vector <ScalarType> &x,
                                                std::vector <ScalarType> &logValue);

/**
 * Computes the exponential map on the hyper-sphere.
 *
 * @param p log-point to exponentiate
 * @param x reference point of the exponential
 * @param expValue exponential value (result)
 */
template <class ScalarType> void sphere_exp_map(const std::vector <ScalarType> &p, const std::vector <ScalarType> &x,
                                                std::vector <ScalarType> &expValue);

/**
 * Computes the weighted average of several points on the hyper-sphere.
 *
 * @param dataPoints List of points on the sphere
 * @param initPoint Starting point for average computation
 * @param weights Weights for each point
 * @param workLogVector Optional work variable (for speed up) to compute the logarithms, to avoid creating it inside this function
 * @param workVector Optional work variable (for speed up), to avoid creating it inside this function
 * @param tol Tolerance for convergence
 */
template <class ScalarType> void ComputeSphericalCentroid(const std::vector < std::vector <ScalarType> > &dataPoints,
                                                          std::vector <ScalarType> &centroidValue, const std::vector <ScalarType> &initPoint,
                                                          const std::vector <ScalarType> &weights, std::vector <ScalarType> *workLogVector = 0,
                                                          std::vector <ScalarType> *workVector = 0, double tol = 1.0e-4);

} // end namespace anima

#include "animaLogExpMapsUnitSphere.hxx"
