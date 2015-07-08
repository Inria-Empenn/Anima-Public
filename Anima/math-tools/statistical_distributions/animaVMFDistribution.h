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
    
} // end of namespace

#include "animaVMFDistribution.hxx"
