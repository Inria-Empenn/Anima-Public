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

    CholeskyDecomposition(const MatrixType &val)
    {
        m_InputMatrix = val;
        m_MatrixSize = val.rows();
        m_LMatrix.set_size(m_MatrixSize, m_MatrixSize);
        m_LMatrix.set_identity();
        m_DMatrix.set_size(m_MatrixSize);
    }

    void PerformDecomposition();
    VectorType SolveLinearSystem(const VectorType &b);
    void Update(const VectorType &x);
    MatrixType Recompose();

private:
    MatrixType m_InputMatrix, m_LMatrix;
    DiagonalType m_DMatrix;
    unsigned int m_MatrixSize;
};

} // end of namespace anima
