#pragma once

#include <vector>
#include "AnimaStatisticalTestsExport.h"

namespace anima
{

// Pair comparison functor
struct pair_decreasing_comparator
{
    bool operator()(const std::pair<unsigned int,double> & f, const std::pair<unsigned int,double> & s)
    {
        return (f.second > s.second);
    }
};

//! In place correction of p-values according to Benjamini Hochberg FDR method
ANIMASTATISTICALTESTS_EXPORT void BHCorrection(std::vector <double> &pvalues);

} // end of namespace anima


