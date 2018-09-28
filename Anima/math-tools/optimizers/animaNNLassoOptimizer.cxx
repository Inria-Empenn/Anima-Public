#include <animaNNLassoOptimizer.h>

namespace anima
{

void NNLassoOptimizer::StartOptimization()
{
    unsigned int numEquations = m_DataMatrix.rows();
    unsigned int numAtoms = m_DataMatrix.cols();

    // Prepare internal variables
    m_BVector.resize(numAtoms);
    for (unsigned int i = 0;i < numAtoms;++i)
    {
        m_BVector[i] = m_L1NormWeight;
        for (unsigned int j = 0;j < numEquations;++j)
            m_BVector[i] -= 2.0 * m_DataMatrix(j,i) * m_Points[j];
    }

    // Now perform successive multiplicative updates
    m_CurrentPosition.SetSize(numAtoms);
    m_CurrentPosition.Fill(1.0);

    ParametersType oldPosition = m_CurrentPosition;
    bool continueLoop = true;
    while (continueLoop)
    {
        double aFactor = 0;
        double cFactor = 0;
        continueLoop = false;

        for (unsigned int i = 0;i < numAtoms;++i)
        {
            for (unsigned int j = 0;j < numAtoms;++j)
            {
                aFactor += m_PositiveSquaredDataMatrix(i,j) * oldPosition[j];
                cFactor += m_NegativeSquaredDataMatrix(i,j) * oldPosition[j];
            }

            m_CurrentPosition[i] = (std::sqrt(m_BVector[i] * m_BVector[i] + 4.0 * aFactor * cFactor) - m_BVector[i]) * oldPosition[i] / (2.0 * aFactor);

            if (!continueLoop)
            {
                if (std::abs(m_CurrentPosition[i] - oldPosition[i]) > 1.0e-6)
                    continueLoop = true;
            }
        }
    }
}

void NNLassoOptimizer::SetDataMatrix(MatrixType &data)
{
    m_DataMatrix = data;
    vnl_matrix <double> squaredDataMatrix = data.transpose() * data;

    m_PositiveSquaredDataMatrix = squaredDataMatrix;
    m_NegativeSquaredDataMatrix = squaredDataMatrix;

    for (unsigned int i = 0;i < squaredDataMatrix.rows();++i)
    {
        for (unsigned int j = 0;j < squaredDataMatrix.cols();++j)
        {
            if (squaredDataMatrix(i,j) < 0)
                m_PositiveSquaredDataMatrix(i,j) = 0;
            else
                m_NegativeSquaredDataMatrix(i,j) = 0;
        }
    }
}

}
