#pragma once
#include "animaDistributionSampling.h"

#include <cmath>
#include <boost/math/distributions/beta.hpp>

#include <animaVectorOperations.h>
#include <animaLogarithmFunctions.h>
#include <animaBaseTensorTools.h>
#include <animaMatrixOperations.h>

#include <itkMacro.h>

namespace anima
{

template <class VectorType>
void SampleFromUniformDistributionOn2Sphere(std::mt19937 &generator, VectorType &resVec)
{
    std::uniform_real_distribution<double> uniDbl(-1.0,1.0);
    double sqSum = 2;
    while (sqSum > 1)
    {
        resVec[0] = uniDbl(generator);
        resVec[1] = uniDbl(generator);

        sqSum = resVec[0] * resVec[0] + resVec[1] * resVec[1];
    }

    double factor = 2.0 * std::sqrt(1.0 - sqSum);
    resVec[0] *= factor;
    resVec[1] *= factor;
    resVec[2] = 2.0 * sqSum - 1.0;
}

template <class T>
unsigned int SampleFromBernoulliDistribution(const T &p, std::mt19937 &generator)
{
    std::bernoulli_distribution bernoulli(p);
    return bernoulli(generator);
}

template <class T>
double SampleFromGaussianDistribution(const T &mean, const T &std, std::mt19937 &generator)
{
    std::normal_distribution<T> normalDist(mean,std);
    return normalDist(generator);
}

template <class VectorType, class ScalarType>
void SampleFromMultivariateGaussianDistribution(const VectorType &mean, const vnl_matrix <ScalarType> &mat, VectorType &resVec,
                                                std::mt19937 &generator, bool isMatCovariance)
{
    unsigned int vectorSize = mat.rows();

    vnl_matrix <ScalarType> stdMatrix = mat;

    if (isMatCovariance)
    {
        vnl_matrix <ScalarType> eVecs(vectorSize,vectorSize);
        vnl_diag_matrix <ScalarType> eVals(vectorSize);

        itk::SymmetricEigenAnalysis < vnl_matrix <ScalarType>, vnl_diag_matrix<ScalarType>, vnl_matrix <ScalarType> > eigenComputer(vectorSize);
        eigenComputer.ComputeEigenValuesAndVectors(mat, eVals, eVecs);

        for (unsigned int i = 0;i < vectorSize;++i)
            eVals[i] = std::sqrt(eVals[i]);

        anima::RecomposeTensor(eVals,eVecs,stdMatrix);
    }

    VectorType tmpVec (resVec);
    for (unsigned int i = 0;i < vectorSize;++i)
        tmpVec[i] = SampleFromGaussianDistribution(0.0,1.0,generator);

    for (unsigned int i = 0;i < vectorSize;++i)
    {
        resVec[i] = mean[i];
        for (unsigned int j = 0;j < vectorSize;++j)
            resVec[i] += stdMatrix(i,j) * tmpVec[j];
    }
}

} // end of namespace anima
