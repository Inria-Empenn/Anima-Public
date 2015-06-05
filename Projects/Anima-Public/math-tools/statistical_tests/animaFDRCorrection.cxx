#include "animaFDRCorrection.h"
#include <algorithm>
#include <iostream>

namespace anima
{
void BHCorrection(std::vector <double> &pvalues)
{
    std::vector < std::pair <unsigned int, double> > indexedPvalues;
    unsigned int numData = pvalues.size();
    for (unsigned int i = 0;i < numData;++i)
        indexedPvalues.push_back(std::pair <unsigned int, double> (i,pvalues[i]));

    std::sort(indexedPvalues.begin(),indexedPvalues.end(),pair_decreasing_comparator());

    // Update p-values
    for (unsigned int i = 0;i < numData;++i)
    {
        double newPvalue = std::min(1.0,indexedPvalues[i].second * numData / (numData - i));
        if (i > 0)
        {
            if (newPvalue > indexedPvalues[i-1].second)
                newPvalue = indexedPvalues[i-1].second;
        }

        indexedPvalues[i].second = newPvalue;
    }

    // Now put back values to their places
    for (unsigned int i = 0;i < numData;++i)
        pvalues[indexedPvalues[i].first] = indexedPvalues[i].second;
}

} // end of namespace anima
