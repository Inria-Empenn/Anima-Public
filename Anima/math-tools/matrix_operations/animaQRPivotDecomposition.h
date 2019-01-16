#pragma once
#include <vnl_matrix.h>
#include <vector>

namespace anima
{

template <typename ScalarType> void QRPivotDecomposition(vnl_matrix <ScalarType> &aMatrix, std::vector <unsigned int> &pivotVector, unsigned int &rank);

} // end namespace anima

#include "animaQRPivotDecomposition.hxx"
