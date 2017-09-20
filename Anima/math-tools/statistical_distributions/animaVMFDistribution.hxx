#pragma once
#include <cmath>

#include "animaVMFDistribution.h"
#include <animaVectorOperations.h>
#include <animaLogarithmFunctions.h>

#include <vnl/vnl_matrix.h>

#include <itkExceptionObject.h>

namespace anima
{
    
template <class VectorType, class ScalarType> double ComputeVMFPdf(const VectorType &v, const VectorType &meanDirection, const ScalarType &kappa)
{
    if (std::abs(anima::ComputeNorm(meanDirection) - 1.0) > 1.0e-6 || std::abs(anima::ComputeNorm(v) - 1.0) > 1.0e-6)
    {
        std::cout << "Norm of mean direction: " << anima::ComputeNorm(meanDirection) << std::endl;
        std::cout << "Norm of evaluated direction: " << anima::ComputeNorm(v) << std::endl;
        throw itk::ExceptionObject(__FILE__, __LINE__,"Von Mises & Fisher distribution requires unitary vectors.",ITK_LOCATION);
    }

    double resVal;

    if (kappa < 1e-4)
    {
        resVal = std::exp(kappa * anima::ComputeScalarProduct(meanDirection, v)) / (4.0 * M_PI);
    }
    else
    {
        double tmpVal = anima::ComputeScalarProduct(meanDirection, v) - 1.0;
        tmpVal *= kappa;

        resVal = std::exp(tmpVal);
        resVal *= kappa;
        tmpVal = 1.0 - std::exp(-2.0 * kappa);
        resVal /= (2.0 * M_PI * tmpVal);
    }

    return resVal;
}

template <class ScalarType, unsigned int Dimension>
double
VMFDistance(const itk::Point <ScalarType, Dimension> &muFirst, const double &kappaFirst,
            const itk::Point <ScalarType, Dimension> &muSec, const double &kappaSec)
{
    double dist = 0;

    double scalarProduct = 0;
    for (unsigned int i = 0;i < Dimension;++i)
        scalarProduct += muFirst[i] * muSec[i];

    if (scalarProduct < -1.0)
        scalarProduct = -1.0;

    if (scalarProduct > 1.0)
        scalarProduct = 1.0;

    double acosValue = std::acos(scalarProduct);

    dist = anima::safe_log(std::max(1.0e-16,kappaSec) / std::max(1.0e-16,kappaFirst));
    dist *= dist;
    dist += acosValue * acosValue;

    return sqrt(dist);
}

template <class ScalarType>
double
GetVonMisesConcentrationMLE(const ScalarType rbar)
{
    double kappa = 0.0;
    
    if (rbar < 0.44639)
    {
        double rbarSq = rbar * rbar;
        kappa = rbar * (12.0 + 6.0 * rbarSq + 5.0 * rbarSq * rbarSq) / 6.0;
    }
    else if (rbar < 0.48070)
        kappa = 0.1 * (rbar - 0.44639) / (0.48070 - 0.44639) + 1.0;
    else if (rbar < 0.51278)
        kappa = 0.1 * (rbar - 0.48070) / (0.51278 - 0.48070) + 1.1;
    else if (rbar < 0.54267)
        kappa = 0.1 * (rbar - 0.51278) / (0.54267 - 0.51278) + 1.2;
    else if (rbar < 0.57042)
        kappa = 0.1 * (rbar - 0.54267) / (0.57042 - 0.54267) + 1.3;
    else if (rbar < 0.59613)
        kappa = 0.1 * (rbar - 0.57042) / (0.59613 - 0.57042) + 1.4;
    else if (rbar < 0.61990)
        kappa = 0.1 * (rbar - 0.59613) / (0.61990 - 0.59613) + 1.5;
    else if (rbar < 0.64183)
        kappa = 0.1 * (rbar - 0.61990) / (0.64183 - 0.61990) + 1.6;
    else if (rbar < 0.66204)
        kappa = 0.1 * (rbar - 0.64183) / (0.66204 - 0.64183) + 1.7;
    else if (rbar < 0.68065)
        kappa = 0.1 * (rbar - 0.66204) / (0.68065 - 0.66204) + 1.8;
    else if (rbar < 0.69777)
        kappa = 0.1 * (rbar - 0.68065) / (0.69777 - 0.68065) + 1.9;
    else if (rbar < 0.71353)
        kappa = 0.1 * (rbar - 0.69777) / (0.71353 - 0.69777) + 2.0;
    else if (rbar < 0.72803)
        kappa = 0.1 * (rbar - 0.71353) / (0.72803 - 0.71353) + 2.1;
    else if (rbar < 0.74138)
        kappa = 0.1 * (rbar - 0.72803) / (0.74138 - 0.72803) + 2.2;
    else if (rbar < 0.75367)
        kappa = 0.1 * (rbar - 0.74138) / (0.75367 - 0.74138) + 2.3;
    else if (rbar < 0.76500)
        kappa = 0.1 * (rbar - 0.75367) / (0.76500 - 0.75367) + 2.4;
    else if (rbar < 0.77545)
        kappa = 0.1 * (rbar - 0.76500) / (0.77545 - 0.76500) + 2.5;
    else if (rbar < 0.78511)
        kappa = 0.1 * (rbar - 0.77545) / (0.78511 - 0.77545) + 2.6;
    else if (rbar < 0.79404)
        kappa = 0.1 * (rbar - 0.78511) / (0.79404 - 0.78511) + 2.7;
    else if (rbar < 0.80231)
        kappa = 0.1 * (rbar - 0.79404) / (0.80231 - 0.79404) + 2.8;
    else
    {
        double irbar = 1.0 - rbar;
        double irbarSq = irbar * irbar;
        kappa = 1.0 / (2.0 * irbar - irbarSq - irbar * irbarSq);
    }
    
    return kappa;
}
    
} // end of namespace anima
