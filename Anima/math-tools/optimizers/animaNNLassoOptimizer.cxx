#include <animaNNLassoOptimizer.h>
#include <animaNNLSOptimizer.h>

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

        double secondTerm = 0.0;
        for (unsigned int j = 0;j < numEquations;++j)
            secondTerm += m_DataMatrix(j,i) * m_Points[j];
        m_BVector[i] -= 2.0 * secondTerm;
    }

    m_CurrentPosition.SetSize(numAtoms);
    m_CurrentPosition.Fill(1.0);

    // First initialize solution from NNLS
    anima::NNLSOptimizer::Pointer nnlsInitializer = anima::NNLSOptimizer::New();
    nnlsInitializer->SetDataMatrix(m_DataMatrix);
    nnlsInitializer->SetPoints(m_Points);
    nnlsInitializer->StartOptimization();

    m_CurrentPosition = nnlsInitializer->GetCurrentPosition();

    double sumWeights = 0.0;
    for (unsigned int i = 0;i < numAtoms;++i)
        sumWeights += m_CurrentPosition[i];

    std::cout << sumWeights << " " << sumWeights / numAtoms << std::endl;
    m_CurrentPosition.Fill(sumWeights / numAtoms);

    // Now perform successive multiplicative updates
    ParametersType oldPosition = m_CurrentPosition;
    bool continueLoop = true;
    double aFactor = 0.0;
    double cFactor = 0.0;
    unsigned int numIter = 0;

    while (continueLoop)
    {
        ++numIter;
        if (numIter % 2000 == 0)
            std::cout << numIter << " " << m_CurrentPosition << std::endl;
        continueLoop = false;

        for (unsigned int i = 0;i < numAtoms;++i)
        {
            aFactor = 0.0;
            cFactor = 0.0;
            for (unsigned int j = 0;j < numAtoms;++j)
            {
                aFactor += m_PositiveSquaredDataMatrix(i,j) * oldPosition[j];
                cFactor += m_NegativeSquaredDataMatrix(i,j) * oldPosition[j];
            }

            m_CurrentPosition[i] = (std::sqrt(m_BVector[i] * m_BVector[i] + 4.0 * aFactor * cFactor) - m_BVector[i]) * oldPosition[i] / (2.0 * aFactor);

            if (!continueLoop)
            {
                double absDiff = std::abs(m_CurrentPosition[i] - oldPosition[i]);
                if (absDiff < 1.0e-6)
                    continue;

                if (absDiff > 1.0e-3 * oldPosition[i])
                    continueLoop = true;
            }
        }

        if (continueLoop)
            oldPosition = m_CurrentPosition;
    }

    std::cout << numIter << std::endl;
}

void NNLassoOptimizer::SetDataMatrix(MatrixType &data)
{
    m_DataMatrix = data;
    vnl_matrix <double> squaredDataMatrix = data.transpose() * data * 2.0;

    m_PositiveSquaredDataMatrix = squaredDataMatrix;
    m_NegativeSquaredDataMatrix = squaredDataMatrix;

    for (unsigned int i = 0;i < squaredDataMatrix.rows();++i)
    {
        for (unsigned int j = 0;j < squaredDataMatrix.cols();++j)
        {
            if (squaredDataMatrix(i,j) < 0)
            {
                m_PositiveSquaredDataMatrix(i,j) = 0;
                m_NegativeSquaredDataMatrix(i,j) = std::abs(m_NegativeSquaredDataMatrix(i,j));
            }
            else
                m_NegativeSquaredDataMatrix(i,j) = 0;
        }
    }
}

}
