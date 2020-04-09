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

}
