#include <animaRootFindingAlgorithms.h>

namespace anima
{

bool CheckRootTolerance::operator()(const double &a, const double &b)
{
    return std::abs(b - a) < m_RootRelativeTolerance * (a + b) / 2.0;
}

double BisectionRootFindingAlgorithm::Optimize()
{
    bool continueLoop = true;
    unsigned int nbIterations = 0;
    
    while (continueLoop)
    {
        ++nbIterations;
        
        double rootValue = (m_LowerBound + m_UpperBound) / 2.0;
        double zeroValue = m_RootFindingFunction(rootValue);
        
        if (m_UseZeroTolerance)
            continueLoop = (std::abs(zeroValue) >= m_ZeroRelativeTolerance);
        
        if (zeroValue < 0.0)
            m_UpperBound = p[0];
        else
            m_LowerBound = p[0];
        
        if (nbIterations >= m_MaximumNumberOfIterations || std::abs(m_UpperBound - m_LowerBound) < xTolRel * (m_LowerBound + m_UpperBound) / 2.0)
            continueLoop = false;
    }
}

}
