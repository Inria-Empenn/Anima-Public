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
    m_CurrentPosition.Fill(0);

    std::vector <unsigned short> treatedIndexes(parametersSize,0);

    std::vector <double> wVector(parametersSize,0);
    std::vector <double> tmpVector(numEquations,0);
    VectorType sPVector(parametersSize);
    sPVector.fill(0);
    VectorType dataPointsP(numEquations);
    dataPointsP.fill(0);
    MatrixType dataMatrixP(numEquations,parametersSize);

    unsigned int numProcessedIndexes = 0;
    std::vector <unsigned int> processedIndexes;

    for (unsigned int i = 0;i < parametersSize;++i)
    {
        for (unsigned int j = 0;j < numEquations;++j)
            wVector[i] += m_DataMatrix(j,i) * m_Points[j];
    }

    bool continueMainLoop = true;
    double previousMaxW = -1;
    while (continueMainLoop)
    {
        double maxW = 0;
        int maxIndex = -1;
        for (unsigned int i = 0;i < parametersSize;++i)
        {
            if (treatedIndexes[i] == 0)
            {
                if (maxIndex < 0)
                {
                    maxW = wVector[i];
                    maxIndex = i;
                }
                else if (maxW < wVector[i])
                {
                    maxW = wVector[i];
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
        treatedIndexes[maxIndex] = 1;
        numProcessedIndexes = this->UpdateProcessedIndexes(treatedIndexes,processedIndexes);
        this->ComputeSPVector(dataMatrixP,dataPointsP,processedIndexes,
                              sPVector,numProcessedIndexes);

        double minSP = sPVector[0];
        for (unsigned int i = 1;i < numProcessedIndexes;++i)
        {
            if (sPVector[i] < minSP)
                minSP = sPVector[i];
        }

        // Starting inner loop
        while (minSP <= 0)
        {
            bool foundOne = false;
            double alpha = 0;
            double indexSelected = 0;
            for (unsigned int i = 0;i < numProcessedIndexes;++i)
            {
                if (sPVector[i] > 0)
                    continue;

                double tmpAlpha = m_CurrentPosition[processedIndexes[i]] / (m_CurrentPosition[processedIndexes[i]] - sPVector[i] + m_EpsilonValue);
                if ((tmpAlpha < alpha)||(!foundOne))
                {
                    alpha = tmpAlpha;
                    indexSelected = processedIndexes[i];
                    foundOne = true;
                }
            }

            for (unsigned int i = 0;i < parametersSize;++i)
                m_CurrentPosition[i] -= alpha * m_CurrentPosition[i];

            // Update positions
            for (unsigned int i = 0;i < numProcessedIndexes;++i)
            {
                m_CurrentPosition[processedIndexes[i]] += alpha * sPVector[i];

                if (std::abs(m_CurrentPosition[processedIndexes[i]]) <= m_EpsilonValue)
                    treatedIndexes[processedIndexes[i]] = 0;
            }

            m_CurrentPosition[indexSelected] = 0;
            treatedIndexes[indexSelected] = 0;

            // Update processed indexes
            numProcessedIndexes = this->UpdateProcessedIndexes(treatedIndexes,processedIndexes);
            if (numProcessedIndexes == 0)
                break;

            this->ComputeSPVector(dataMatrixP,dataPointsP,processedIndexes,
                                  sPVector,numProcessedIndexes);

            minSP = sPVector[0];
            for (unsigned int i = 1;i < numProcessedIndexes;++i)
            {
                if (sPVector[i] < minSP)
                    minSP = sPVector[i];
            }
        }

        m_CurrentPosition.Fill(0);
        for (unsigned int i = 0;i < numProcessedIndexes;++i)
            m_CurrentPosition[processedIndexes[i]] = sPVector[i];

        if (numProcessedIndexes == parametersSize)
        {
            continueMainLoop = false;
            continue;
        }

        for (unsigned int i = 0;i < numEquations;++i)
        {
            tmpVector[i] = m_Points[i];
            for (unsigned int j = 0;j < parametersSize;++j)
                tmpVector[i] -= m_DataMatrix(i,j) * m_CurrentPosition[j];
        }

        for (unsigned int i = 0;i < parametersSize;++i)
        {
            wVector[i] = 0;
            for (unsigned int j = 0;j < numEquations;++j)
                wVector[i] += m_DataMatrix(j,i) * tmpVector[j];
        }
    }
}

unsigned int NNLSOptimizer::UpdateProcessedIndexes(std::vector <unsigned short> &treatedIndexes,
                                                   std::vector <unsigned int> &processedIndexes)
{
    unsigned int numProcessedIndexes = 0;
    unsigned int parametersSize = m_DataMatrix.cols();

    for (unsigned int i = 0;i < parametersSize;++i)
    {
        if (treatedIndexes[i] > 0)
            numProcessedIndexes++;
    }

    processedIndexes.resize(numProcessedIndexes);
    unsigned int pos = 0;
    for (unsigned int i = 0;i < parametersSize;++i)
    {
        if (treatedIndexes[i] > 0)
        {
            processedIndexes[pos] = i;
            ++pos;
        }
    }

    return numProcessedIndexes;
}

void NNLSOptimizer::ComputeSPVector(MatrixType &dataMatrixP, VectorType &dataPointsP,
                                    std::vector<unsigned int> &processedIndexes,
                                    VectorType &sPVector,unsigned int numProcessedIndexes)
{
    unsigned int numEquations = m_DataMatrix.rows();

    dataMatrixP.set_size(numEquations,numProcessedIndexes);
    dataPointsP.set_size(numEquations);

    for (unsigned int i = 0;i < numEquations;++i)
    {
        for (unsigned int j = 0;j < numProcessedIndexes;++j)
            dataMatrixP(i,j) = m_DataMatrix(i,processedIndexes[j]);

        dataPointsP[i] = m_Points[i];
    }

    sPVector = vnl_qr<double>(dataMatrixP).solve(dataPointsP);
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
