#include <animaBVLSOptimizer.h>
#include <vnl/algo/vnl_qr.h>

namespace anima
{

void BVLSOptimizer::StartOptimization()
{
//    std::cout << "Optim with data " << m_DataMatrix << std::endl;
//    std::cout << "And vector " << m_Points << std::endl;
//    std::cout << "With lower bounds " << m_LowerBounds << std::endl;
//    std::cout << "With upper bounds " << m_UpperBounds << std::endl;

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
//        std::cout << "Finding T*..." << std::endl;
        for (unsigned int i = 0;i < parametersSize;++i)
        {
            double testVal = m_WVector[i] * m_ParametersAtBoundsVector[i];
//            std::cout << i << " " << m_WVector[i] << " " << m_ParametersAtBoundsVector[i] << " " << testVal << std::endl;
            if (testVal > maxVal)
            {
                maxVal = testVal;
                maxIndex = i;
            }
        }

//        std::cout << "T* value " << maxIndex << " " << maxVal << std::endl;

        // Set it to F
        m_ParametersAtBoundsVector[maxIndex] = 0;

        bool continueInternalLoop = true;
        while (continueInternalLoop)
        {
            // Compute b' and A'
            unsigned int numReducedParameters = 0;
            for (unsigned int i = 0;i < parametersSize;++i)
                numReducedParameters += (m_ParametersAtBoundsVector[i] == 0);

//            std::cout << "Number of reduced parameters " << numReducedParameters << std::endl;
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
                        m_ReducedDataMatrix(j,pos) = m_DataMatrix(j,k);
                        ++pos;
                    }
                }
            }

//            std::cout << "Reduced data " << m_ReducedDataMatrix << std::endl;
//            std::cout << "b' vector " << m_bPrimeVector << std::endl;

            // Compute reduced solution
            m_ReducedSolution = vnl_qr <double> (m_ReducedDataMatrix).solve(m_bPrimeVector);

//            std::cout << m_ReducedSolution << std::endl;

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

//                std::cout << "Reduced ok, done " << m_CurrentPosition << std::endl;
            }
            else
            {
                // If not ok, compute alphas and update sets
                // Init minAlpha to 2.0, will be replaced at first out bound
                // since it should be between 0.0 and 1.0
                double minAlpha = 2.0;
                unsigned int minIndex = 0;

//                std::cout << "Reduced not ok " << std::endl;
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

//                std::cout << "Min alpha and index " << minAlpha << " " << minIndex << " " << maxIndex << std::endl;

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
//                std::cout << "Bounding again some paramters" << std::endl;
                for (unsigned int i = 0;i < parametersSize;++i)
                {
                    if (m_ParametersAtBoundsVector[i] != 0)
                        continue;

//                    std::cout << i << " " << m_CurrentPosition[i] << " " << minAlpha << " " << m_ReducedSolution[pos];

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

//                    std::cout << " " << m_CurrentPosition[i] << " " << m_LowerBounds[i] << " " << m_UpperBounds[i];
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

//                    std::cout << " " << m_ParametersAtBoundsVector[i] << std::endl;

                    ++pos;
                }
            }

            wRequired = true;
        }
    }

//    std::cout << "Done optimizing" << std::endl;
}

void BVLSOptimizer::InitializeSolutionByProjection()
{
    m_CurrentPosition = vnl_qr <double> (m_DataMatrix).solve(m_Points);

    unsigned int parametersSize = m_CurrentPosition.size();
    m_ParametersAtBoundsVector.resize(parametersSize);

    for (unsigned int i = 0;i < parametersSize;++i)
    {
        m_ParametersAtBoundsVector[i] = 0;
        if (m_CurrentPosition[i] <= m_LowerBounds[i])
        {
            m_CurrentPosition[i] = m_LowerBounds[i];
            m_ParametersAtBoundsVector[i] = 1;
        }
        else if (m_CurrentPosition[i] >= m_UpperBounds[i])
        {
            m_CurrentPosition[i] = m_UpperBounds[i];
            m_ParametersAtBoundsVector[i] = -1;
        }
        else
        {
            // Warm start as in Stark and Parker
            m_CurrentPosition[i] = (m_LowerBounds[i] + m_UpperBounds[i]) / 2.0;
        }
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
            m_TmpVector[i] -= m_DataMatrix(i,j) * m_CurrentPosition[j];
    }

    for (unsigned int i = 0;i < numParameters;++i)
    {            
        m_WVector[i] = 0.0;

        for (unsigned int j = 0;j < numEqs;++j)
            m_WVector[i] += m_DataMatrix(j,i) * m_TmpVector[j];
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

    if ((freeSetFull) || (encounteredL && allNegativeWsForL) || (encounteredU && allPositiveWsForU))
        return false;

    return true;
}

double BVLSOptimizer::GetCurrentResidual()
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
