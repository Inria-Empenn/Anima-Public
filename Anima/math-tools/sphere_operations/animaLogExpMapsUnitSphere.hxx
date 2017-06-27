#pragma once
#include "animaLogExpMapsUnitSphere.h"
#include <cmath>

namespace anima
{

template <class ScalarType> void sphere_log_map(const std::vector <ScalarType> &p, const std::vector <ScalarType> &x,
                                                std::vector <ScalarType> &logValue)
{
    unsigned int sizeVec = p.size();
    if (logValue.size() != sizeVec)
        logValue.resize(sizeVec);

    for (unsigned int i = 0;i < sizeVec;++i)
        logValue[i] = 0;

    double dotProd = 0;
    for (unsigned int i = 0;i < sizeVec;++i)
        dotProd += p[i]*x[i];

    if (dotProd >= 1)
        return;

    if (dotProd <= -1 + 1.0e-8)
        dotProd = -1 + 1.0e-8;

    for (unsigned int i = 0;i < sizeVec;++i)
        logValue[i] = (p[i] - x[i]*dotProd)*std::acos(dotProd)/std::sqrt(1.0 - dotProd*dotProd);
}

template <class ScalarType> void sphere_exp_map(const std::vector <ScalarType> &p, const std::vector <ScalarType> &x,
                                                std::vector <ScalarType> &expValue)
{
    unsigned int sizeVec = p.size();
    if (expValue.size() != sizeVec)
        expValue.resize(sizeVec);

    double normP = 0;
    for (unsigned int i = 0;i < sizeVec;++i)
        normP += p[i]*p[i];

    normP = std::sqrt(normP);

    if (normP == 0)
    {
        expValue = x;
        return;
    }

    for (unsigned int i = 0;i < sizeVec;++i)
        expValue[i] = std::cos(normP)*x[i] + p[i]*std::sin(normP)/normP;
}

template <class ScalarType> void ComputeSphericalCentroid(const std::vector < std::vector <ScalarType> > &dataPoints,
                                                          std::vector <ScalarType> &centroidValue, const std::vector <ScalarType> &initPoint,
                                                          const std::vector <ScalarType> &weights, std::vector <ScalarType> *workLogVector,
                                                          std::vector <ScalarType> *workVector, double tol)
{
    unsigned int sizeData = dataPoints.size();
    unsigned int sizeVec = initPoint.size();
    centroidValue = initPoint;

    bool internalLogVector = false;
    bool internalVector = false;
    if (!workLogVector)
    {
        workLogVector = new std::vector <ScalarType> (sizeVec);
        internalLogVector = true;
    }

    if (workLogVector->size() != sizeVec)
        workLogVector->resize(sizeVec);

    if (!workVector)
    {
        workVector = new std::vector <ScalarType> (sizeVec);
        internalVector = true;
    }

    if (workVector->size() != sizeVec)
        workVector->resize(sizeVec);

    double sumWeights = 0;
    for (unsigned int i = 0;i < sizeData;++i)
        sumWeights += weights[i];

    double diffNewOld = tol + 1;

    while (diffNewOld > tol)
    {
        for (unsigned int j = 0;j < sizeVec;++j)
            (*workLogVector)[j] = 0;

        for (unsigned int i = 0;i < sizeData;++i)
        {
            sphere_log_map(dataPoints[i],centroidValue,(*workVector));

            for (unsigned int j = 0;j < sizeVec;++j)
                (*workLogVector)[j] += weights[i]*(*workVector)[j];
        }

        for (unsigned int j = 0;j < sizeVec;++j)
            (*workLogVector)[j] /= sumWeights;

        sphere_exp_map((*workLogVector),centroidValue,(*workVector));

        diffNewOld = 0;

        for (unsigned int j = 0;j < sizeVec;++j)
            diffNewOld += (centroidValue[j]-(*workVector)[j])*(centroidValue[j]-(*workVector)[j]);

        diffNewOld = std::sqrt(diffNewOld);

        centroidValue = (*workVector);
    }

    if (internalLogVector)
        delete workLogVector;

    if (internalVector)
        delete workVector;
}

} //end namespace anima
