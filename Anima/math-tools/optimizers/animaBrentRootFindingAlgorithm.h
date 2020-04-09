#pragma once
#include "AnimaOptimizersExport.h"

#include <animaBaseRootFindingAlgorithm.h>

namespace anima
{

class ANIMAOPTIMIZERS_EXPORT BrentRootFindingAlgorithm : public BaseRootFindingAlgorithm
{
public:
    using Superclass = BaseRootFindingAlgorithm;
    using ParametersType = Superclass::ParametersType;
    double Optimize() ITK_OVERRIDE;
};

} // end namespace anima
