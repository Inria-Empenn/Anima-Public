#include <animaMatrixDecompositions.h>

namespace anima
{

void CholeskyDecomposition::PerformDecomposition()
{
    for (unsigned int j = 0;j < m_MatrixSize;++j)
    {
        m_DMatrix[j] = m_InputMatrix(j,j);

        if (j > 0)
            for (unsigned int k = 0;k < j;++k)
                m_DMatrix[j] -= m_LMatrix(j,k) * m_LMatrix(j,k) * m_DMatrix[k];

        for (unsigned int i = j + 1;i < m_MatrixSize;++i)
        {
            double tmpVal = m_InputMatrix(i,j);

            if (j > 0)
                for (unsigned int k = 0;k < j;++k)
                    tmpVal -= m_LMatrix(i,k) * m_LMatrix(j,k) * m_DMatrix[k];

            m_LMatrix(i,j) = tmpVal / m_DMatrix[j];
        }
    }
}

void CholeskyDecomposition::SolveLinearSystem(VectorType &b)
{
    for (unsigned int i = 0;i < m_MatrixSize;++i)
    {
        double intermediateValue = b[i];

        if (i > 0)
            for (unsigned int j = 0;j < i;++j)
                intermediateValue -= m_LMatrix(i,j) * b[j];

        b[i] = intermediateValue;
    }

    for (int i = m_MatrixSize - 1;i >= 0;--i)
    {
        double intermediateValue = b[i] / m_DMatrix[i];

        if (i < m_MatrixSize - 1)
            for (unsigned int j = i + 1;j < m_MatrixSize;++j)
                intermediateValue -= m_LMatrix(j,i) * b[j];

        b[i] = intermediateValue;
    }
}

void CholeskyDecomposition::Update(VectorType &x)
{
    double a = 1.0;

    for (unsigned int i = 0;i < m_MatrixSize;++i)
    {
        double p = x[i];
        double oldDValue = m_DMatrix[i];
        m_DMatrix[i] += a * p * p;
        double b = p * a / m_DMatrix[i];
        a *= oldDValue / m_DMatrix[i];

        if (i < m_MatrixSize - 1)
        {
            for (unsigned int j = i + 1;j < m_MatrixSize;++j)
            {
                x[j] -= p * m_LMatrix(j,i);
                m_LMatrix(j,i) += b * x[j];
            }
        }
    }
}

void CholeskyDecomposition::Recompose()
{
    for (unsigned int i = 0;i < m_MatrixSize;++i)
    {
    	for (unsigned int j = i;j < m_MatrixSize;++j)
        {
        	double tmpVal = 0.0;
            for (unsigned int k = 0;k <= i;++k)
                tmpVal += m_LMatrix(i,k) * m_DMatrix[k] * m_LMatrix(j,k);

            m_InputMatrix(i,j) = tmpVal;

            if (i != j)
            	m_InputMatrix(j,i) = tmpVal;
        }
    }
}

}
