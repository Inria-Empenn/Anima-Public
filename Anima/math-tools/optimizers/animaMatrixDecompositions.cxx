#include <animaMatrixDecompositions.h>

namespace anima
{

void CholeskyDecomposition::PerformDecomposition()
{
    for (unsigned int j = 0;j < m_MatrixSize;++j)
    {
        m_DMatrix[j] = m_InputMatrix(j,j);

        if (j > 0)
            for (unsigned int k = 0;k <= j - 1;++k)
                m_DMatrix[j] -= m_LMatrix(j,k) * m_LMatrix(j,k) * m_DMatrix[k];

        for (unsigned int i = j + 1;i < m_MatrixSize;++i)
        {
            double tmpVal = m_InputMatrix(i,j);

            if (j > 0)
                for (unsigned int k = 0;k <= j - 1;++k)
                    tmpVal -= m_LMatrix(i,k) * m_LMatrix(j,k) * m_DMatrix[k];

            m_LMatrix(i,j) = tmpVal / m_DMatrix[j];
        }
    }
}

CholeskyDecomposition::VectorType CholeskyDecomposition::SolveLinearSystem(const VectorType &b)
{
    VectorType resVal(m_MatrixSize);

    for (unsigned int i = 0;i < m_MatrixSize;++i)
    {
        double intermediateValue = b[i];

        if (i > 0)
            for (unsigned int j = 0;j <= i - 1;++j)
                intermediateValue -= m_LMatrix(i,j) * resVal[j];

        resVal[i] = intermediateValue;
    }

    for (int i = m_MatrixSize - 1;i >= 0;--i)
    {
        double intermediateValue = resVal[i] / m_DMatrix[i];

        if (i < m_MatrixSize - 1)
            for (unsigned int j = i + 1;j < m_MatrixSize;++j)
                intermediateValue -= m_LMatrix(j,i) * resVal[j];

        resVal[i] = intermediateValue;
    }

    return resVal;
}

void CholeskyDecomposition::Update(const VectorType &x)
{
    double a = 1.0;
    VectorType w = x;

    for (unsigned int i = 0;i < m_MatrixSize;++i)
    {
        double p = w[i];
        double d0 = m_DMatrix[i];
        m_DMatrix[i] += a * p * p;
        double b = p * a / m_DMatrix[i];
        a *= d0 / m_DMatrix[i];

        for (unsigned int j = i + 1;j < m_MatrixSize;++j)
        {
            w[j] -= p * m_LMatrix(j,i);
            m_LMatrix(j,i) += b * w[j];
        }
    }
}

}
