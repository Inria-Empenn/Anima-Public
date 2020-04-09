#include "animaBaseRootFindingAlgorithm.h"
#include <itkMacro.h>

namespace anima
{

bool CheckRootTolerance::operator()(const double &a, const double &b)
{
    return std::abs(b - a) < m_RootRelativeTolerance * (a + b) / 2.0;
}

double RootFindingFunctionBoostBridge::operator()(const double &x)
{
    m_ParametersVector[0] = x;
    return m_RootFindingFunction->GetValue(m_ParametersVector);
}

void BaseRootFindingAlgorithm::SetFunctionValueAtInitialLowerBound(const double &val)
{
    m_FunctionValueAtInitialLowerBound = val;
    m_ProvidedFunctionValueAtInitialLowerBound = true;
}

void BaseRootFindingAlgorithm::SetFunctionValueAtInitialUpperBound(const double &val)
{
    m_FunctionValueAtInitialUpperBound = val;
    m_ProvidedFunctionValueAtInitialUpperBound = true;
}

}
