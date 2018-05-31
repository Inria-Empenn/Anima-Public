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

//! Maximum likelihood estimation of the concentration parameter of the 2D von Mises distribution according to Mardia, Statistics of Directional Data, 1972.
template <class ScalarType>
double GetVonMisesConcentrationMLE(const ScalarType rbar);
    
} // end of namespace

#include "animaVMFDistribution.hxx"
