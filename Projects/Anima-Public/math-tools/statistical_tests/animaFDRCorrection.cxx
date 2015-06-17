#include "animaFDRCorrection.h"
#include <algorithm>
#include <iostream>

namespace anima
{

void BHCorrection(std::vector <double> &pvalues, double qValue)
{
    std::vector < std::pair <unsigned int, double> > indexedPvalues;
    unsigned int numData = pvalues.size();
    for (unsigned int i = 0;i < numData;++i)
        indexedPvalues.push_back(std::pair <unsigned int, double> (i,pvalues[i]));

    std::sort(indexedPvalues.begin(),indexedPvalues.end(),pair_increasing_comparator());

    // Update p-values
    unsigned int breakingIndex = 0;
    for (unsigned int i = 0;i < numData;++i)
    {
        if (indexedPvalues[i].second > qValue * (i + 1) / numData)
        {
            breakingIndex = i;
            break;
        }
    }

    // Now put back values to their places
    for (unsigned int i = 0;i < numData;++i)
        pvalues[indexedPvalues[i].first] = (i < breakingIndex);
}

} // end of namespace anima
