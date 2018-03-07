#pragma once

#include "animaChiDistribution.h"
#include <animaARBBesselFunctions.h>
#include <animaARBKummerFunctions.h>

namespace anima
{
    
inline double factorial(unsigned int x)
{
    return (x <= 1) ? 1 : x * factorial(x - 1);
}

inline double double_factorial(unsigned int x)
{
    return (x <= 1) ? 1 : x * double_factorial(x - 2);
}

double LowerBound(const unsigned int N, const double k0ValueSq, const double betaNSq)
{
    double xiValue = 2.0 * N - betaNSq * k0ValueSq;
    double insideValue = 2.0 * N / xiValue - 1.0;
    return std::sqrt(std::max(insideValue, 0.0));
}

double XiFunction(const double thetaSq, const unsigned int N, const double k1Value, const double betaNSq)
{
    double resVal = 2.0 * N + thetaSq - betaNSq * k1Value * k1Value;
    
    if (resVal > 1.0)
    {
        std::cerr << "The xi function cannot take values greater than 1." << std::endl;
        std::cout << thetaSq << " " << k1Value << " " << anima::GetKummerM(-thetaSq / 2.0, -0.5, N) << std::endl;
        exit(-1);
    }
    
    return resVal;
}

double GFunction(const double thetaSq, const double rSq, const unsigned int N, const double k1Value, const double betaNSq)
{
    double insideValue = anima::XiFunction(thetaSq, N, k1Value, betaNSq) * (1.0 + rSq) - 2 * N;
    return std::sqrt(std::max(insideValue, 0.0));
}

double KFunction(const double theta, const double r, const unsigned int N, const double betaNSq, double &k1Value)
{
    double thetaSq = theta * theta;
    double rSq = r * r;
    k1Value = anima::GetKummerM(-thetaSq / 2.0, -0.5, N);
    double k2Value = anima::GetKummerM(-thetaSq / 2.0, 0.5, N + 1);
    double gValue = anima::GFunction(thetaSq, rSq, N, k1Value, betaNSq);
    double num = gValue * (gValue - theta);
    double denom = theta * (1.0 + rSq) * (1.0 - betaNSq / (2.0 * N) * k1Value * k2Value) - gValue;
    return theta - num / denom;
}

double FixedPointFinder(const double r, const unsigned int N, double &k1Value, const unsigned int maximumNumberOfIterations, const double epsilon)
{
    double k0Value = anima::GetKummerM(0.0, -0.5, N);
    double k0ValueSq = k0Value * k0Value;
    double betaN = std::sqrt(M_PI / 2.0) * anima::double_factorial(2 * N - 1) / (std::pow(2, N - 1) * anima::factorial(N - 1));
    double betaNSq = betaN * betaN;
    
    double lb = anima::LowerBound(N, k0ValueSq, betaNSq);
    if (r <= lb)
        return 0.0;
    
    unsigned int counter = maximumNumberOfIterations;
    
    double t0 = r - lb;
    double t1 = anima::KFunction(t0, r, N, betaNSq, k1Value);

    while (std::abs(t0 - t1) > epsilon)
    {
        t0 = t1;
        t1 = anima::KFunction(t0, r, N, betaNSq, k1Value);
        --counter;

        if (counter == 0)
            break;
    }
    
    return t1;
}

double GFunction2(const double etaSq, const double mSq, const double sigmaSq, const unsigned int N, const double k1Value, const double betaNSq)
{
    double insideValue = mSq + (anima::XiFunction(etaSq / sigmaSq, N, k1Value, betaNSq) - 2 * N) * sigmaSq;
    return std::sqrt(std::max(insideValue, 0.0));
}

double KFunction2(const double eta, const double m, const double sigma, const unsigned int N, const double betaNSq, double &k1Value)
{
    double mSq = m * m;
    double etaSq = eta * eta;
    double sigmaSq = sigma * sigma;
    double thetaSq = etaSq / sigmaSq;
    k1Value = anima::GetKummerM(-thetaSq / 2.0, -0.5, N);
    double k2Value = anima::GetKummerM(-thetaSq / 2.0, 0.5, N + 1);
    double gValue = anima::GFunction2(etaSq, mSq, sigmaSq, N, k1Value, betaNSq);
    double num = gValue * (gValue - eta);
    double denom = eta * (1.0 - betaNSq / (2.0 * N) * k1Value * k2Value) - gValue;
    return eta - num / denom;
}

double FixedPointFinder2(const double m, const double sigma, const unsigned int N, double &k1Value, const unsigned int maximumNumberOfIterations, const double epsilon)
{
    unsigned int counter = maximumNumberOfIterations;
    double betaN = std::sqrt(M_PI / 2.0) * anima::double_factorial(2 * N - 1) / (std::pow(2, N - 1) * anima::factorial(N - 1));
    double betaNSq = betaN * betaN;
    double delta = betaN * sigma - m;
    
    if (delta == 0)
        return 0;
    
    double mCorrected = (delta > 0) ? betaN * sigma + delta : m;
    
    double t0 = mCorrected;
    double t1 = anima::KFunction2(t0, mCorrected, sigma, N, betaNSq, k1Value);
    
    while (std::abs(t0 - t1) > epsilon)
    {
        t0 = t1;
        t1 = anima::KFunction2(t0, mCorrected, sigma, N, betaNSq, k1Value);
        --counter;
        
        if (counter == 0)
            break;
    }
    
    return (delta > 0) ? -t1 : t1;
}

void GetRiceParameters(const std::vector<double> &samples, const std::vector<double> &weights, double &location, double &scale)
{
    unsigned int sampleSize = samples.size();
    double k1Value = 0;
    double meanValue = 0;
    double sumWeights = 0, sumSqWeights = 0;
    
    for (unsigned int i = 0;i < sampleSize;++i)
    {
        double weightValue = weights[i];
        meanValue += weightValue * samples[i];
        sumWeights += weightValue;
        sumSqWeights += weightValue * weightValue;
    }
    
    meanValue /= sumWeights;
    
    if (scale == 0)
    {
        double sigmaValue = 0;
        
        for (unsigned int i = 0;i < sampleSize;++i)
            sigmaValue += weights[i] * (samples[i] - meanValue) * (samples[i] - meanValue);
        
        sigmaValue /= (sumWeights * sumWeights - sumSqWeights);
        sigmaValue *= sumWeights;
        sigmaValue = std::sqrt(sigmaValue);
        
        double rValue = meanValue / sigmaValue;
        double thetaValue = FixedPointFinder(rValue, 1, k1Value);
        k1Value = anima::GetKummerM(-thetaValue * thetaValue / 2.0, -0.5, 1);
        
        scale = sigmaValue / std::sqrt(anima::XiFunction(thetaValue * thetaValue, 1, k1Value, M_PI / 2.0));
        location = thetaValue * scale;
        return;
    }
    
    location = FixedPointFinder2(meanValue, scale, 1, k1Value);
}

double GetRiceCDF(const double x, const double location, const double scale)
{
    double a = location / scale;
    double b = x / scale;
    return 1.0 - anima::GetMarcumQ(1, a, b);
}
    
} // end of namespace anima
