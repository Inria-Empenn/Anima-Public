#include <animaNNLSOptimizer.h>
#include <vnl/algo/vnl_qr.h>

namespace anima
{

const double NNLSOptimizer::m_EpsilonValue = 1.0e-8;

void NNLSOptimizer::StartOptimization()
{
    unsigned int parametersSize = m_DataMatrix.cols();
    unsigned int numEquations = m_DataMatrix.rows();

    if ((numEquations != m_Points.GetSize())||(numEquations == 0)||(parametersSize == 0))
        itkExceptionMacro("Wrongly sized inputs to NNLS, aborting");

    m_CurrentPosition.SetSize(parametersSize);
    m_TreatedIndexes.resize(parametersSize);
    m_ProcessedIndexes.clear();
    m_WVector.resize(parametersSize);

    unsigned int numProcessedIndexes = 0;

    for (unsigned int i = 0;i < parametersSize;++i)
    {
        m_TreatedIndexes[i] = 0;
        m_CurrentPosition[i] = 0.0;
        m_WVector[i] = 0;
        for (unsigned int j = 0;j < numEquations;++j)
            m_WVector[i] += m_DataMatrix(j,i) * m_Points[j];
    }

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

        std::fill(m_WVector.begin(),m_WVector.end(),0.0);
        for (unsigned int i = 0;i < numEquations;++i)
        {
            double tmpValue = m_Points[i];
            for (unsigned int j = 0;j < parametersSize;++j)
                tmpValue -= m_DataMatrix(i,j) * m_CurrentPosition[j];

            for (unsigned int j = 0;j < parametersSize;++j)
                m_WVector[j] += m_DataMatrix(i,j) * tmpValue;
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

    m_DataMatrixP.set_size(numEquations,numProcessedIndexes);
    m_DataPointsP.set_size(numEquations);

    for (unsigned int i = 0;i < numEquations;++i)
    {
        for (unsigned int j = 0;j < numProcessedIndexes;++j)
            m_DataMatrixP(i,j) = m_DataMatrix(i,m_ProcessedIndexes[j]);

        m_DataPointsP[i] = m_Points[i];
    }

    m_SPVector = vnl_qr <double> (m_DataMatrixP).solve(m_DataPointsP);
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
