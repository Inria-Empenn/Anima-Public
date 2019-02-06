#pragma once
#include <vnl/vnl_matrix.h>
#include <vector>

namespace anima
{

/**
 * Computes QR decomposition with pivoting following Golub et al. Matrix compputations. 1996 (algorithm 5.4.1)
 * Computes aMatrix = Q * R * tranposePivot
 * Everything gets stored back into aMatrix and housBetaValues (R is the upper trianuglar part of aMatrix)
 */
template <typename ScalarType> void QRPivotDecomposition(vnl_matrix <ScalarType> &aMatrix, std::vector <unsigned int> &transposePivotVector,
                                                         std::vector <ScalarType> &houseBetaValues, unsigned int &rank);

/**
 * From the results of QR decomposition, compute QtB from the obtained matrix in QRPivotDecomposition, beta values and provided vector B
 */
template <typename ScalarType> void GetQtBFromQRDecomposition(vnl_matrix <ScalarType> &qrMatrix, vnl_vector <ScalarType> &bVector,
                                                              std::vector <ScalarType> &houseBetaValues, unsigned int rank);

} // end namespace anima

#include "animaQRPivotDecomposition.hxx"
