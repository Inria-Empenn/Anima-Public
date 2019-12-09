#pragma once
#include <vnl/vnl_matrix.h>
#include <vector>

namespace anima
{

/**
 * Computes QR decomposition following Golub et al. Matrix computations. 1996 (algorithm 5.2.2)
 * Computes aMatrix = Q * R
 * Only the R matrix gets stored back into aMatrix. If bVector is not empty, it is replaced by Q^t bVector
 */
template <typename ScalarType> void QRGivensDecomposition(vnl_matrix <ScalarType> &aMatrix, vnl_vector <ScalarType> &bVector);

/**
 * Computes QR decomposition with pivoting following Golub et al. Matrix computations. 1996 (algorithm 5.4.1)
 * Computes aMatrix * pivot = Q * R
 * Everything gets stored back into aMatrix and housBetaValues (R is the upper triangular part of aMatrix)
 * Pivot is stored as pivotVector allowing to construct pivot by setting pivot(pivotVector(i),i) = 1, the rest being null
 */
template <typename ScalarType> void QRPivotDecomposition(vnl_matrix <ScalarType> &aMatrix, std::vector <unsigned int> &pivotVector,
                                                         std::vector <ScalarType> &houseBetaValues, unsigned int &rank);

/**
 * From the results of QR decomposition, compute QtB from the obtained matrix in QRPivotDecomposition, beta values and provided vector B
 */
template <typename ScalarType> void GetQtBFromQRPivotDecomposition(vnl_matrix <ScalarType> &qrMatrix, vnl_vector <ScalarType> &bVector,
                                                                   std::vector <ScalarType> &houseBetaValues, unsigned int rank);

/**
 * From the results of QR decomposition, compute Q from the obtained matrix in QRPivotDecomposition, beta values
 */
template <typename ScalarType> void GetQMatrixQRPivoDecomposition(vnl_matrix <ScalarType> &qrMatrix, std::vector <ScalarType> &houseBetaValues,
                                                                  vnl_matrix <ScalarType> &qMatrix, unsigned int rank);

} // end namespace anima

#include "animaQRPivotDecomposition.hxx"
