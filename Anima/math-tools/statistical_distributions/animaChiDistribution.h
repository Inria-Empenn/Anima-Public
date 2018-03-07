#pragma once

#include <vector>

namespace anima
{

//! Implementation of Koay and Basser, Analytically exact correction scheme for signal extraction from noisy magnitude MR signals, Journal of Magnetic Resonance (2006).
double LowerBound(const unsigned int N, const double k0ValueSq, const double betaNSq);
double XiFunction(const double thetaSq, const unsigned int N, const double k1Value, const double betaNSq);
double GFunction(const double thetaSq, const double rSq, const unsigned int N, const double k1Value, const double betaNSq);
double KFunction(const double theta, const double r, const unsigned int N, const double betaNSq, double &k1Value);
double FixedPointFinder(const double r, const unsigned int N, double &k1Value, const unsigned int maximumNumberOfIterations = 500, const double epsilon = 1.0e-9);

//! Implementation of Koay, Ozarslan and Basser, A signal transformational framework for breaking the noise floor and its applications in MRI, Journal of Magnetic Resonance (2009).
double GFunction2(const double etaSq, const double mSq, const double sigmaSq, const unsigned int N, const double k1Value, const double betaNSq);
double KFunction2(const double eta, const double m, const double sigma, const unsigned int N, const double betaNSq, double &k1Value);
double FixedPointFinder2(const double m, const double sigma, const unsigned int N, double &k1Value, const unsigned int maximumNumberOfIterations = 500, const double epsilon = 1.0e-9);

//! Implementation of both Koay technique for Rice parameter estimation
void GetRiceParameters(const std::vector<double> &samples, const std::vector<double> &weights, double &location, double &scale);

//! In-house function for evaluating the cumulative distribution function of a Rice distribution based on Gaussian quadrature integral approximation of the Marcum Q function.
double GetRiceCDF(const double x, const double location, const double scale);
    
} // end of namespace

#include "animaChiDistribution.hxx"
