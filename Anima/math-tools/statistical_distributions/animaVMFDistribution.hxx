#pragma once
#include <cmath>

#include "animaVMFDistribution.h"
#include <animaVectorOperations.h>
#include <animaLogarithmFunctions.h>

#include <boost/math/distributions/beta.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
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
    
} // end of namespace anima
