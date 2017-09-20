#pragma once

#include <vector>
#include <itkPoint.h>

namespace anima
{
    
template <class VectorType, class ScalarType>
double ComputeVMFPdf(const VectorType &v, const VectorType &meanDirection, const ScalarType &kappa);

template <class ScalarType, unsigned int Dimension>
double VMFDistance(const itk::Point <ScalarType,Dimension> &muFirst, const double &kappaFirst,
                   const itk::Point <ScalarType,Dimension> &muSec, const double &kappaSec);

//! Implementation of an approximation of the MLE of the concentration parameter of the von Mises Fisher distribution according to Mardia, Directional statistics, 1972.
template <class ScalarType>
double GetVonMisesConcentrationMLE(const ScalarType rbar);
    
} // end of namespace

#include "animaVMFDistribution.hxx"
