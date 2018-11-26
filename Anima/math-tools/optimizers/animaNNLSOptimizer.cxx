#include <animaNNLSOptimizer.h>
#include <vnl_qr.h>

namespace anima
{

const double NNLSOptimizer::m_EpsilonValue = 1.0e-12;

void NNLSOptimizer::StartOptimization()
{
    unsigned int parametersSize = m_DataMatrix.cols();
    unsigned int numEquations = m_DataMatrix.rows();

    if ((numEquations != m_Points.size())||(numEquations == 0)||(parametersSize == 0))
        itkExceptionMacro("Wrongly sized inputs to NNLS, aborting");

    m_CurrentPosition.SetSize(parametersSize);
    m_CurrentPosition.Fill(0.0);
    m_TreatedIndexes.resize(parametersSize);
    std::fill(m_TreatedIndexes.begin(), m_TreatedIndexes.end(),0);
    m_ProcessedIndexes.clear();
    m_WVector.resize(parametersSize);

    unsigned int numProcessedIndexes = 0;

    this->ComputeWVector();

    bool continueMainLoop = true;
    double previousMaxW = -1;
    while (continueMainLoop)
    {
        double maxW = 0;
        int maxIndex = -1;
        for (unsigned int i = 0;i < parametersSize;++i)
        {
            if (m_TreatedIndexes[i] == 0)
            {
                if (maxIndex < 0)
                {
                    maxW = m_WVector[i];
                    maxIndex = i;
                }
                else if (maxW < m_WVector[i])
                {
                    maxW = m_WVector[i];
                    maxIndex = i;
                }
            }
        }

        if (maxW <= m_EpsilonValue)
        {
            continueMainLoop = false;
            continue;
        }

        previousMaxW = maxW;
        m_TreatedIndexes[maxIndex] = 1;
        m_ProcessedIndexes.push_back(maxIndex);
        numProcessedIndexes = m_ProcessedIndexes.size();
        this->ComputeSPVector();

        double minSP = m_SPVector[0];
        for (unsigned int i = 1;i < numProcessedIndexes;++i)
        {
            if (m_SPVector[i] < minSP)
                minSP = m_SPVector[i];
        }

        // Starting inner loop
        while (minSP <= 0)
        {
            bool foundOne = false;
            double alpha = 0;
            double indexSelected = 0;
            for (unsigned int i = 0;i < numProcessedIndexes;++i)
            {
                if (m_SPVector[i] > 0)
                    continue;

                double tmpAlpha = m_CurrentPosition[m_ProcessedIndexes[i]] / (m_CurrentPosition[m_ProcessedIndexes[i]] - m_SPVector[i] + m_EpsilonValue);
                if ((tmpAlpha < alpha)||(!foundOne))
                {
                    alpha = tmpAlpha;
                    indexSelected = m_ProcessedIndexes[i];
                    foundOne = true;
                }
            }

            for (unsigned int i = 0;i < parametersSize;++i)
                m_CurrentPosition[i] -= alpha * m_CurrentPosition[i];

            // Update positions
            for (unsigned int i = 0;i < numProcessedIndexes;++i)
            {
                m_CurrentPosition[m_ProcessedIndexes[i]] += alpha * m_SPVector[i];

                if (std::abs(m_CurrentPosition[m_ProcessedIndexes[i]]) <= m_EpsilonValue)
                    m_TreatedIndexes[m_ProcessedIndexes[i]] = 0;
            }

            m_CurrentPosition[indexSelected] = 0;
            m_TreatedIndexes[indexSelected] = 0;

            // Update processed indexes
            numProcessedIndexes = this->UpdateProcessedIndexes();
            if (numProcessedIndexes == 0)
                break;

            this->ComputeSPVector();

            minSP = m_SPVector[0];
            for (unsigned int i = 1;i < numProcessedIndexes;++i)
            {
                if (m_SPVector[i] < minSP)
                    minSP = m_SPVector[i];
            }
        }

        m_CurrentPosition.Fill(0);
        for (unsigned int i = 0;i < numProcessedIndexes;++i)
            m_CurrentPosition[m_ProcessedIndexes[i]] = m_SPVector[i];

        if (numProcessedIndexes == parametersSize)
        {
            continueMainLoop = false;
            continue;
        }

        this->ComputeWVector();
    }
}

void NNLSOptimizer::ComputeWVector()
{
    unsigned int parametersSize = m_DataMatrix.cols();
    unsigned int numEquations = m_DataMatrix.rows();

    m_WVector.resize(parametersSize);

    std::fill(m_WVector.begin(),m_WVector.end(),0.0);
    if (!m_SquaredProblem)
    {
        for (unsigned int i = 0;i < numEquations;++i)
        {
            double tmpValue = m_Points[i];
            for (unsigned int j = 0;j < parametersSize;++j)
                tmpValue -= m_DataMatrix(i,j) * m_CurrentPosition[j];

            for (unsigned int j = 0;j < parametersSize;++j)
                m_WVector[j] += m_DataMatrix(i,j) * tmpValue;
        }
    }
    else
    {
        for (unsigned int i = 0;i < numEquations;++i)
        {
            m_WVector[i] = m_Points[i];
            for (unsigned int j = 0;j < parametersSize;++j)
                m_WVector[i] -= m_DataMatrix(i,j) * m_CurrentPosition[j];
        }
    }
}

unsigned int NNLSOptimizer::UpdateProcessedIndexes()
{
    unsigned int numProcessedIndexes = 0;
    unsigned int parametersSize = m_DataMatrix.cols();

    for (unsigned int i = 0;i < parametersSize;++i)
        numProcessedIndexes += m_TreatedIndexes[i];

    m_ProcessedIndexes.resize(numProcessedIndexes);
    unsigned int pos = 0;
    for (unsigned int i = 0;i < parametersSize;++i)
    {
        if (m_TreatedIndexes[i] != 0)
        {
            m_ProcessedIndexes[pos] = i;
            ++pos;
        }
    }

    return numProcessedIndexes;
}

void NNLSOptimizer::ComputeSPVector()
{
    unsigned int numEquations = m_DataMatrix.rows();
    unsigned int numProcessedIndexes = m_ProcessedIndexes.size();

    if (!m_SquaredProblem)
    {
        m_DataMatrixP.set_size(numEquations,numProcessedIndexes);
        m_SPVector.set_size(numProcessedIndexes);
        m_SPVector.fill(0.0);

        // We do not square on the fly since we are not allowed to, squared matrix has a chance of being too bad
        for (unsigned int i = 0;i < numEquations;++i)
        {
            for (unsigned int j = 0;j < numProcessedIndexes;++j)
                m_DataMatrixP(i,j) = m_DataMatrix(i,m_ProcessedIndexes[j]);
        }

        m_SPVector = vnl_qr <double> (m_DataMatrixP).solve(m_Points);
    }
    else
    {
        m_DataMatrixP.set_size(numProcessedIndexes,numProcessedIndexes);
        m_SPVector.set_size(numProcessedIndexes);
        m_DataMatrixP.fill(0.0);
        m_SPVector.fill(0.0);

        for (unsigned int i = 0;i < numProcessedIndexes;++i)
        {
            unsigned int iIndex = m_ProcessedIndexes[i];
            m_SPVector[i] = m_Points[iIndex];
            for (unsigned int j = i;j < numProcessedIndexes;++j)
            {
                unsigned int jIndex = m_ProcessedIndexes[j];
                m_DataMatrixP(i,j) = m_DataMatrix(iIndex,jIndex);

                if (i != j)
                    m_DataMatrixP(j,i) = m_DataMatrixP(i,j);
            }
        }

        m_CholeskySolver.SetInputMatrix(m_DataMatrixP);
        m_CholeskySolver.PerformDecomposition();
        m_CholeskySolver.SolveLinearSystemInPlace(m_SPVector);
    }
}

double NNLSOptimizer::GetCurrentResidual()
{
    double residualValue = 0;

    for (unsigned int i = 0;i < m_DataMatrix.rows();++i)
    {
        double tmpVal = 0;
        for (unsigned int j = 0;j < m_DataMatrix.cols();++j)
            tmpVal += m_DataMatrix(i,j) * m_CurrentPosition[j];

        residualValue += (tmpVal - m_Points[i]) * (tmpVal - m_Points[i]);
    }

    return residualValue;
}

}
