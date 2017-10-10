#include "AnimaOptimizersExport.h"

namespace anima
{

// Utility functions to bound and rescale values, used in all compartments
namespace levenberg
{
ANIMAOPTIMIZERS_EXPORT double UnboundValue(double x, double lowerBound, double upperBound);
ANIMAOPTIMIZERS_EXPORT double ComputeBoundedValue(double x, double &inputSign, double lowerBound, double upperBound);
ANIMAOPTIMIZERS_EXPORT double BoundedDerivativeAddOn(double x, double inputSign, double lowerBound, double upperBound);
}

} // end namespace anima
