#pragma once

#include <vnl_matrix.h>
#include <vnl_diag_matrix.h>
#include <vnl_vector.h>

#include "AnimaOptimizersExport.h"

namespace anima
{

/**
 * @brief Cholesky decomposition: decomposes a symmetric matrix A in the form L D L^T,
 * where L is lower triangular and D diagonal. May be used to solve efficiently any linear system
 * using the solve methods. Refer to Gill, Golub et al. Methods for modifying matrix factorizations. 1974
 */
class ANIMAOPTIMIZERS_EXPORT CholeskyDecomposition
{
public:
    typedef vnl_matrix<double> MatrixType;
    typedef vnl_diag_matrix<double> DiagonalType;
    typedef vnl_vector<double> VectorType;

    CholeskyDecomposition() {}

    CholeskyDecomposition(const unsigned int matrixDimension, const double epsilon)
    {
        m_InputMatrix.set_size(matrixDimension, matrixDimension);
        m_InputMatrix.fill(0.0);
        m_InputMatrix.fill_diagonal(epsilon);
        m_MatrixSize = matrixDimension;
        m_LMatrix.set_size(matrixDimension, matrixDimension);
        m_LMatrix.set_identity();
        m_DMatrix.set_size(matrixDimension);
        m_DMatrix.fill_diagonal(epsilon);
    }

    CholeskyDecomposition(const MatrixType &val)
    {
        this->SetInputMatrix(val);
    }

    //! Performs LDL decomposition from input matrix
    void PerformDecomposition();

    //! Recompose input matrix from LDL decomposition
    void Recompose();

    //! Solves linear system Ax=b and outputs result in new variable
    VectorType &SolveLinearSystem(const VectorType &b);

    //! Solves linear system Ax=b and outputs result in input b variable
    void SolveLinearSystemInPlace(VectorType &b);

    //! Update decomposition with x so that LDL matches A + x x^T
    void Update(const VectorType &x);

    //! Set input matrix to decompose
    void SetInputMatrix(const MatrixType &matrix);

    //! Get condition number from decomposition
    double GetConditionNumber();

    //! Get input matrix, useful if decomposition updated through rank one method
    MatrixType &GetInputMatrix() {return m_InputMatrix;}

    void SetLMatrix(const MatrixType &val) {m_LMatrix = val;}
    MatrixType &GetLMatrix() {return m_LMatrix;}
    void SetDMatrix(const DiagonalType &val) {m_DMatrix = val;}
    DiagonalType &GetDMatrix() {return m_DMatrix;}

private:
    MatrixType m_InputMatrix, m_LMatrix;
    DiagonalType m_DMatrix;
    unsigned int m_MatrixSize;

    VectorType m_ProblemSolution;
    VectorType m_WorkVector;
};

} // end of namespace anima
