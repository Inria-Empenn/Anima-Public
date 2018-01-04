#include <animaMatchingPursuitOptimizer.h>

namespace anima
{

void MatchingPursuitOptimizer::StartOptimization()
{
    unsigned int dictionarySize = m_DataMatrix.cols();
    unsigned int numEquations = m_DataMatrix.rows();

    if (dictionarySize < m_MaximalNumberOfWeights)
        m_MaximalNumberOfWeights = dictionarySize;

    m_CurrentPosition.SetSize(dictionarySize);
    m_CurrentPosition.Fill(0);

    std::vector <unsigned int> indexesTaken;
    unsigned int numIndexesDone = 0;

    m_CurrentResiduals.resize(numEquations);
    for (unsigned int i = 0;i < numEquations;++i)
        m_CurrentResiduals[i] = m_Points[i];

    std::vector <double> dictionaryNorms(dictionarySize,0.0);
    for (unsigned int i = 0;i < dictionarySize;++i)
    {
        for (unsigned int j = 0;j < numEquations;++j)
            dictionaryNorms[i] += m_DataMatrix(j,i) * m_DataMatrix(j,i);

        dictionaryNorms[i] = std::sqrt(dictionaryNorms[i]);
    }

    while(numIndexesDone < m_MaximalNumberOfWeights)
    {
        double maxDotProductDictionary = 0;
        int maxIndex = -1;
        for (unsigned int i = 0;i < dictionarySize;++i)
        {
            if (!this->CheckIndex(indexesTaken,i))
                continue;

            double dotProductDictionary = 0;
            for (unsigned int j = 0;j < numEquations;++j)
                dotProductDictionary += m_DataMatrix(j,i) * m_CurrentResiduals[j];

            dotProductDictionary /= dictionaryNorms[i];
            if (dotProductDictionary > maxDotProductDictionary)
            {
                maxDotProductDictionary = dotProductDictionary;
                maxIndex = i;
            }
        }

        if (maxIndex < 0)
            break;

        if (maxIndex >= m_IgnoredIndexesUpperBound)
            ++numIndexesDone;

        indexesTaken.push_back((unsigned int)maxIndex);
        m_CurrentPosition[maxIndex] = maxDotProductDictionary / dictionaryNorms[maxIndex];
        for (unsigned int i = 0;i < numEquations;++i)
            m_CurrentResiduals[i] -= m_CurrentPosition[maxIndex] * m_DataMatrix(i,maxIndex);
    }
}

double MatchingPursuitOptimizer::GetCurrentResidual()
{
    unsigned int numEquations = m_DataMatrix.rows();

    double sqResidual = 0;
    for (unsigned int i = 0;i < numEquations;++i)
        sqResidual += m_CurrentResiduals[i] * m_CurrentResiduals[i];

    return sqResidual;
}

bool MatchingPursuitOptimizer::CheckIndex(std::vector <unsigned int> &indexesTaken, unsigned int index)
{
    for (unsigned int i = 0;i < indexesTaken.size();++i)
    {
        if (index == indexesTaken[i])
            return false;
    }

    return true;
}

}
