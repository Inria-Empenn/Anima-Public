#pragma once
#include "AnimaOptimizersExport.h"

#include <animaBaseRootFindingAlgorithm.h>

namespace anima
{

class ANIMAOPTIMIZERS_EXPORT TOMS748RootFindingAlgorithm : public BaseRootFindingAlgorithm
{
public:
    using Superclass = BaseRootFindingAlgorithm;
    using ParametersType = Superclass::ParametersType;
    double Optimize() ITK_OVERRIDE;
};

} // end namespace anima
