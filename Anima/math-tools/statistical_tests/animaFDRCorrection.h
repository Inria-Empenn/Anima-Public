#pragma once

#include <vector>
#include "AnimaStatisticalTestsExport.h"

namespace anima
{
    
    // Pair comparison functor
    struct pair_increasing_comparator
    {
        bool operator()(const std::pair<unsigned int,double> & f, const std::pair<unsigned int,double> & s)
        {
            return (f.second < s.second);
        }
    };
    
    /**
     * In place correction of p-values according to Benjamini Hochberg FDR method
     * Output is a thresholded list at the specified q-value
     * Eq. (1) of Y. Benjamini and Y. Hochberg. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing.
     * Journal of the Royal Statistical Society. Series B (Methodological)
     * Vol. 57, No. 1 (1995), pp. 289-300
     */
    ANIMASTATISTICALTESTS_EXPORT void BHCorrection(std::vector <double> &pvalues, double qValue);
    ANIMASTATISTICALTESTS_EXPORT void BYCorrection(std::vector <double> &pvalues, double qValue);
    
} // end of namespace anima


