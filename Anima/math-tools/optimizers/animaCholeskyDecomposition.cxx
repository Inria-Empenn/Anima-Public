#include <animaCholeskyDecomposition.h>
#include <iostream>
#include <limits>

namespace anima
{

void CholeskyDecomposition::SetInputMatrix(const MatrixType &val)
{
    m_InputMatrix = val;
    m_MatrixSize = val.rows();
}

void CholeskyDecomposition::PerformDecomposition()
{
    m_LMatrix.set_size(m_MatrixSize, m_MatrixSize);
    m_LMatrix.set_identity();
    m_DMatrix.set_size(m_MatrixSize);
\
    for (unsigned int i = 0;i < m_MatrixSize;++i)
    {
        m_DMatrix[i] = m_InputMatrix(i,i);
        for (unsigned int j = 0;j < i;++j)
            m_DMatrix[i] -= m_LMatrix(i,j) * m_LMatrix(i,j) * m_DMatrix[j];

        for (unsigned int j = i + 1;j < m_MatrixSize;++j)
        {
            double diffVal = 0;
            for (unsigned int k = 0;k < i;++k)
                diffVal += m_LMatrix(j,k) * m_LMatrix(i,k) * m_DMatrix[k];

            m_LMatrix(j,i) = (m_InputMatrix(j,i) - diffVal) / m_DMatrix[i];
        }
    }
}

double CholeskyDecomposition::GetConditionNumber()
{
    double minVal = m_DMatrix[0];
    double maxVal = m_DMatrix[0];

    for (unsigned int i = 1;i < m_MatrixSize;++i)
    {
        double absVal = std::abs(m_DMatrix[i]);
        if (minVal > absVal)
            minVal = absVal;

        if (maxVal < absVal)
            maxVal = absVal;
    }

    if (minVal == 0.0)
        return std::numeric_limits<double>::max();

    return maxVal / minVal;
}

CholeskyDecomposition::VectorType &CholeskyDecomposition::SolveLinearSystem(const VectorType &b)
{
    m_ProblemSolution = b;
    SolveLinearSystemInPlace(m_ProblemSolution);
    return m_ProblemSolution;
}

void CholeskyDecomposition::SolveLinearSystemInPlace(VectorType &b)
{
    for (unsigned int i = 1;i < m_MatrixSize;++i)
    {
        for (unsigned int j = 0;j < i;++j)
            b[i] -= m_LMatrix(i,j) * b[j];
    }

    for (unsigned int i = 0;i < m_MatrixSize;++i)
        b[i] /= m_DMatrix[i];

    for (int i = m_MatrixSize - 2;i >= 0;--i)
    {
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
