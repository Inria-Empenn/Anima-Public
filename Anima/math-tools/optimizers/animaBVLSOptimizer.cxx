#include <animaBVLSOptimizer.h>
#include <vnl/algo/vnl_qr.h>

namespace anima
{

void BVLSOptimizer::StartOptimization()
{
    unsigned int parametersSize = m_DataMatrix.cols();
    unsigned int numEquations = m_DataMatrix.rows();

    if ((numEquations != m_Points.size())||(numEquations == 0)||(parametersSize == 0))
        itkExceptionMacro("Wrongly sized inputs to BVLS, aborting");

    m_bPrimeVector.set_size(numEquations);
    m_CurrentPosition.SetSize(parametersSize);
    m_ParametersAtBoundsVector.resize(parametersSize);

    this->InitializeSolutionByProjection();

    bool continueLoop = true;
    bool wRequired = true;

    while (continueLoop)
    {
        if (wRequired)
            this->ComputeWVector();
        continueLoop = this->TestKuhnTuckerConvergence();
        if (!continueLoop)
            break;

        // Find t*
        double maxVal = 0.0;
        unsigned int maxIndex = 0;
        for (unsigned int i = 0;i < parametersSize;++i)
        {
            double testVal = m_WVector[i] * m_ParametersAtBoundsVector[i];
            if (testVal > maxVal)
            {
                maxVal = testVal;
                maxIndex = i;
            }
        }

        // Set it to F
        m_ParametersAtBoundsVector[maxIndex] = 0;

        bool continueInternalLoop = true;
        while (continueInternalLoop)
        {
            // Compute b' and A'
            unsigned int numReducedParameters = 0;
            for (unsigned int i = 0;i < parametersSize;++i)
                numReducedParameters += (m_ParametersAtBoundsVector[i] == 0);

            if (numReducedParameters == 0)
                break;

            m_ReducedDataMatrix.set_size(numEquations,numReducedParameters);
            for (unsigned int j = 0;j < numEquations;++j)
            {
                m_bPrimeVector[j] = m_Points[j];
                unsigned int pos = 0;
                for (unsigned int k = 0;k < parametersSize;++k)
                {
                    if (m_ParametersAtBoundsVector[k] != 0)
                        m_bPrimeVector[j] -= m_DataMatrix(j,k) * m_CurrentPosition[k];
                    else
                    {
                        m_ReducedDataMatrix[j][pos] = m_DataMatrix[j][k];
                        ++pos;
                    }
                }
            }

            // Compute reduced solution
            m_ReducedSolution = vnl_qr <double> (m_ReducedDataMatrix).solve(m_bPrimeVector);

            // Test reduced solution
            bool reducedOk = true;
            unsigned int pos = 0;
            for (unsigned int i = 0;i < parametersSize;++i)
            {
                if (m_ParametersAtBoundsVector[i] != 0)
                    continue;

                if ((m_ReducedSolution[pos] <= m_LowerBounds[i])||
                        (m_ReducedSolution[pos] >= m_UpperBounds[i]))
                {
                    reducedOk = false;
                    break;
                }

                ++pos;
            }

            if (reducedOk)
            {
                continueInternalLoop = false;
                pos = 0;
                for (unsigned int i = 0;i < parametersSize;++i)
                {
                    if (m_ParametersAtBoundsVector[i] != 0)
                        continue;

                    m_CurrentPosition[i] = m_ReducedSolution[pos];
                    ++pos;
                }
            }
            else
            {
                // If not ok, compute alphas and update sets
                // Init minAlpha to 2.0, will be replaced at first out bound
                // since it should be between 0.0 and 1.0
                double minAlpha = 2.0;
                unsigned int minIndex = 0;

                pos = 0;
                for (unsigned int i = 0;i < parametersSize;++i)
                {
                    if (m_ParametersAtBoundsVector[i] != 0)
                        continue;

                    if ((m_ReducedSolution[pos] <= m_LowerBounds[i])||
                            (m_ReducedSolution[pos] >= m_UpperBounds[i]))
                    {
                        double denom = m_ReducedSolution[pos] - m_CurrentPosition[i];
                        double alpha = 0.0;
                        if (m_ReducedSolution[pos] >= m_UpperBounds[i])
                            alpha = m_UpperBounds[i] - m_CurrentPosition[i];
                        else
                            alpha = m_LowerBounds[i] - m_CurrentPosition[i];

                        alpha /= denom;
                        alpha = std::max(0.0,alpha);

                        if (alpha < minAlpha)
                        {
                            minAlpha = alpha;
                            minIndex = i;
                        }
                    }

                    ++pos;
                }

                if ((minAlpha == 0.0) && (minIndex == maxIndex))
                {
                    m_WVector[maxIndex] = 0.0;
                    wRequired = false;
                    continueInternalLoop = false;

                    if (m_CurrentPosition[maxIndex] == m_LowerBounds[maxIndex])
                        m_ParametersAtBoundsVector[maxIndex] = 1;
                    else if (m_CurrentPosition[maxIndex] == m_UpperBounds[maxIndex])
                        m_ParametersAtBoundsVector[maxIndex] = -1;

                    continue;
                }

                pos = 0;
                for (unsigned int i = 0;i < parametersSize;++i)
                {
                    if (m_ParametersAtBoundsVector[i] != 0)
                        continue;

                    if (i != minIndex)
                    {
                        double addonValue = minAlpha * (m_ReducedSolution[pos] - m_CurrentPosition[i]);
                        m_CurrentPosition[i] += addonValue;
                    }
                    else
                    {
                        if (m_ReducedSolution[pos] <= m_LowerBounds[i])
                            m_CurrentPosition[i] = m_LowerBounds[i];
                        else
                            m_CurrentPosition[i] = m_UpperBounds[i];
                    }

                    if (m_CurrentPosition[i] <= m_LowerBounds[i])
                    {
                        m_ParametersAtBoundsVector[i] = 1;
                        m_CurrentPosition[i] = m_LowerBounds[i];
                    }
                    else if (m_CurrentPosition[i] >= m_UpperBounds[i])
                    {
                        m_ParametersAtBoundsVector[i] = -1;
                        m_CurrentPosition[i] = m_UpperBounds[i];
                    }

                    ++pos;
                }
            }

            wRequired = true;
        }
    }
}

void BVLSOptimizer::InitializeSolutionByProjection()
{
    unsigned int parametersSize = m_DataMatrix.cols();
    m_CurrentPosition.SetSize(parametersSize);
    m_ParametersAtBoundsVector.resize(parametersSize);
    std::fill(m_ParametersAtBoundsVector.begin(), m_ParametersAtBoundsVector.end(),0);

    unsigned int countBounded = 0;
    unsigned int diffCountBounded = 1;
    unsigned int numEquations = m_DataMatrix.rows();

    while ((diffCountBounded != 0) && (countBounded < parametersSize))
    {
        unsigned int numParameters = parametersSize - countBounded;

        m_ReducedDataMatrix.set_size(numEquations,numParameters);
        for (unsigned int j = 0;j < numEquations;++j)
        {
            m_bPrimeVector[j] = m_Points[j];
            unsigned int pos = 0;
            for (unsigned int k = 0;k < parametersSize;++k)
            {
                if (m_ParametersAtBoundsVector[k] != 0)
                    m_bPrimeVector[j] -= m_DataMatrix[j][k] * m_CurrentPosition[k];
                else
                {
                    m_ReducedDataMatrix[j][pos] = m_DataMatrix[j][k];
                    ++pos;
                }
            }
        }

        m_ReducedSolution = vnl_qr <double> (m_ReducedDataMatrix).solve(m_bPrimeVector);
        diffCountBounded = 0;
        unsigned int pos = 0;
        for (unsigned int i = 0;i < parametersSize;++i)
        {
            if (m_ParametersAtBoundsVector[i] != 0)
                continue;

            if (m_ReducedSolution[pos] <= m_LowerBounds[i])
            {
                m_CurrentPosition[i] = m_LowerBounds[i];
                m_ParametersAtBoundsVector[i] = 1;
                ++diffCountBounded;
            }
            else if (m_ReducedSolution[pos] >= m_UpperBounds[i])
            {
                m_CurrentPosition[i] = m_UpperBounds[i];
                m_ParametersAtBoundsVector[i] = -1;
                ++diffCountBounded;
            }
            else
                m_CurrentPosition[i] = m_ReducedSolution[pos];

            ++pos;
        }

        countBounded += diffCountBounded;
    }
}

void BVLSOptimizer::ComputeWVector()
{
    unsigned int numEqs = m_Points.size();
    unsigned int numParameters = m_DataMatrix.cols();
    m_TmpVector.resize(numEqs);
    m_WVector.resize(numParameters);

    for (unsigned int i = 0;i < numEqs;++i)
    {
        m_TmpVector[i] = m_Points[i];
        for (unsigned int j = 0;j < numParameters;++j)
            m_TmpVector[i] -= m_DataMatrix[i][j] * m_CurrentPosition[j];
    }

    for (unsigned int i = 0;i < numParameters;++i)
    {            
        m_WVector[i] = 0.0;

        for (unsigned int j = 0;j < numEqs;++j)
            m_WVector[i] += m_DataMatrix[j][i] * m_TmpVector[j];
    }
}

bool BVLSOptimizer::TestKuhnTuckerConvergence()
{
    bool freeSetFull = true;
    bool allNegativeWsForL = true;
    bool allPositiveWsForU = true;
    bool encounteredL = false;
    bool encounteredU = false;

    unsigned int parametersSize = m_DataMatrix.cols();
    for (unsigned int i = 0;i < parametersSize;++i)
    {
        if (m_ParametersAtBoundsVector[i] == 0)
            continue;

        freeSetFull = false;
        if (m_ParametersAtBoundsVector[i] == 1)
        {
            encounteredL = true;
            if (m_WVector[i] > 0.0)
                allNegativeWsForL = false;
        }
        else
        {
            encounteredU = true;
            if (m_WVector[i] < 0.0)
                allPositiveWsForU = false;
        }
    }

    bool returnValue = !freeSetFull;
    if (returnValue)
    {
        if (encounteredL && encounteredU)
            returnValue = !(allNegativeWsForL && allPositiveWsForU);
        else if (encounteredL)
            returnValue = !allNegativeWsForL;
        else if (encounteredU)
            returnValue = !allPositiveWsForU;
    }

    return returnValue;
}

double BVLSOptimizer::GetCurrentResidual()
{
    double residualValue = 0;

    for (unsigned int i = 0;i < m_DataMatrix.rows();++i)
    {
        double tmpVal = 0;
        for (unsigned int j = 0;j < m_DataMatrix.cols();++j)
            tmpVal += m_DataMatrix[i][j] * m_CurrentPosition[j];

        residualValue += (tmpVal - m_Points[i]) * (tmpVal - m_Points[i]);
    }

    return residualValue;
}

}
