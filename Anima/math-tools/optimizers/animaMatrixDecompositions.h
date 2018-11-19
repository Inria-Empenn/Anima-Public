#pragma once

#include <vnl_matrix.h>
#include <vnl_diag_matrix.h>
#include <vnl_vector.h>

#include "AnimaOptimizersExport.h"

namespace anima
{

class ANIMAOPTIMIZERS_EXPORT CholeskyDecomposition
{
public:
    typedef vnl_matrix<double> MatrixType;
    typedef vnl_diag_matrix<double> DiagonalType;
    typedef vnl_vector<double> VectorType;

    CholeskyDecomposition() {}

    CholeskyDecomposition(const unsigned int n, const double val)
    {
        m_InputMatrix.set_size(n, n);
        m_InputMatrix.fill(0.0);
        m_InputMatrix.fill_diagonal(val);
        m_MatrixSize = n;
        m_LMatrix.set_size(n, n);
        m_LMatrix.set_identity();
        m_DMatrix.set_size(n);
        m_DMatrix.fill_diagonal(val);
    }

    CholeskyDecomposition(const MatrixType &val)
    {
        m_InputMatrix = val;
        m_MatrixSize = val.rows();
        m_LMatrix.set_size(m_MatrixSize, m_MatrixSize);
        m_LMatrix.set_identity();
        m_DMatrix.set_size(m_MatrixSize);
    }

    void PerformDecomposition();
    void Recompose();
    void SolveLinearSystem(VectorType &b);
    void Update(VectorType &x);

    MatrixType & GetInputMatrix() {return m_InputMatrix;};
    void SetLMatrix(const MatrixType &val) {m_LMatrix = val;};
    MatrixType & GetLMatrix() {return m_LMatrix;};
    void SetDMatrix(const DiagonalType &val) {m_DMatrix = val;};
    DiagonalType & GetDMatrix() {return m_DMatrix;};

private:
    MatrixType m_InputMatrix, m_LMatrix;
    DiagonalType m_DMatrix;
    unsigned int m_MatrixSize;
};

} // end of namespace anima
