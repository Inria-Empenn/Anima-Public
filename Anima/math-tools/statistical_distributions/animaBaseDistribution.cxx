#include "animaBaseDistribution.h"

namespace anima
{

double BaseDistribution::GetLogDensity(const double &x)
{
    return std::log(this->GetDensity(x));
}
    
} // end of namespace anima
