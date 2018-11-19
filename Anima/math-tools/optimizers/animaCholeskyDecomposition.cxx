#include <animaCholeskyDecomposition.h>

namespace anima
{

void CholeskyDecomposition::SetInputMatrix(const MatrixType &val)
{
    m_InputMatrix = val;
    m_MatrixSize = val.rows();
    m_LMatrix.set_size(m_MatrixSize, m_MatrixSize);
    m_LMatrix.set_identity();
    m_DMatrix.set_size(m_MatrixSize);
}

void CholeskyDecomposition::PerformDecomposition()
{
    for (unsigned int j = 0;j < m_MatrixSize;++j)
    {
        m_DMatrix[j] = m_InputMatrix(j,j);

        for (unsigned int k = 0;k < j;++k)
            m_DMatrix[j] -= m_LMatrix(j,k) * m_LMatrix(j,k) * m_DMatrix[k];

        for (unsigned int i = j + 1;i < m_MatrixSize;++i)
        {
            double tmpVal = m_InputMatrix(i,j);

            for (unsigned int k = 0;k < j;++k)
                tmpVal -= m_LMatrix(i,k) * m_LMatrix(j,k) * m_DMatrix[k];

            m_LMatrix(i,j) = tmpVal / m_DMatrix[j];
        }
    }
}

CholeskyDecomposition::VectorType &CholeskyDecomposition::SolveLinearSystem(const VectorType &b)
{
    m_ProblemSolution.set_size(m_MatrixSize);
    for (unsigned int i = 0;i < m_MatrixSize;++i)
    {
        m_ProblemSolution[i] = b[i];

        for (unsigned int j = 0;j < i;++j)
            m_ProblemSolution[i] -= m_LMatrix(i,j) * m_ProblemSolution[j];
    }

    for (int i = m_MatrixSize - 1;i >= 0;--i)
    {
        m_ProblemSolution[i] /= m_DMatrix[i];

        for (unsigned int j = i + 1;j < m_MatrixSize;++j)
            m_ProblemSolution[i] -= m_LMatrix(j,i) * m_ProblemSolution[j];
    }

    return m_ProblemSolution;
}

void CholeskyDecomposition::SolveLinearSystemInPlace(VectorType &b)
{
    b.set_size(m_MatrixSize);
    for (unsigned int i = 1;i < m_MatrixSize;++i)
    {
        for (unsigned int j = 0;j < i;++j)
            b[i] -= m_LMatrix(i,j) * b[j];
    }

    for (int i = m_MatrixSize - 1;i >= 0;--i)
    {
        b[i] /= m_DMatrix[i];

        for (unsigned int j = i + 1;j < m_MatrixSize;++j)
            b[i] -= m_LMatrix(j,i) * b[j];
    }
}

void CholeskyDecomposition::Update(const VectorType &x)
{
    double a = 1.0;
    m_WorkVector = x;

    for (unsigned int i = 0;i < m_MatrixSize;++i)
    {
        double p = m_WorkVector[i];
        double oldDValue = m_DMatrix[i];
        m_DMatrix[i] += a * p * p;
        double b = p * a / m_DMatrix[i];
        a *= oldDValue / m_DMatrix[i];

        for (unsigned int j = i + 1;j < m_MatrixSize;++j)
        {
            m_WorkVector[j] -= p * m_LMatrix(j,i);
            m_LMatrix(j,i) += b * m_WorkVector[j];
        }
    }
}

void CholeskyDecomposition::Recompose()
{
    for (unsigned int i = 0;i < m_MatrixSize;++i)
    {
        for (unsigned int j = i;j < m_MatrixSize;++j)
        {
            m_InputMatrix(i,j) = 0.0;
            for (unsigned int k = 0;k <= i;++k)
                m_InputMatrix(i,j) += m_LMatrix(i,k) * m_DMatrix[k] * m_LMatrix(j,k);

            if (i != j)
                m_InputMatrix(j,i) = m_InputMatrix(i,j);
        }
    }
}

}
