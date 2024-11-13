#pragma once
#include "AnimaMathsOptimizersExport.h"

#include <animaBaseRootFindingAlgorithm.h>

namespace anima
{

class ANIMAMATHSOPTIMIZERS_EXPORT DekkerRootFindingAlgorithm : public BaseRootFindingAlgorithm
{
public:
    using Superclass = BaseRootFindingAlgorithm;
    using ParametersType = Superclass::ParametersType;
    double Optimize() ITK_OVERRIDE;
};

} // end namespace anima
