#include <animaBaseMCMCost.h>
#include <animaMCMConstants.h>

namespace anima
{

BaseMCMCost::BaseMCMCost()
{
    m_B0Value = 1;
    m_SigmaSquare = 1;

    m_SmallDelta = anima::DiffusionSmallDelta;
    m_BigDelta = anima::DiffusionBigDelta;
}

} // end namespace anima
