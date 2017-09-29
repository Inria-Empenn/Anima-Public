#pragma once
#include "animaSphereOperations.h"
#include <cmath>

namespace anima
{

template <class ScalarType> void GetSphereEvenSampling(std::vector < std::vector <ScalarType> > &spherePoints,
                                                       unsigned int numSamples)
{
    double numLatitudes = std::sqrt(numSamples * M_PI / 8.0);
    double oldNumLatitudes = numLatitudes - 1;

    while (std::abs(numLatitudes - oldNumLatitudes) > 1.0e-8)
    {
        oldNumLatitudes = numLatitudes;
        numLatitudes = numSamples * std::sin(M_PI / (4.0 * numLatitudes)) / 2.0;
    }

    unsigned int roundedNumLatitudes = std::round(numLatitudes);
    spherePoints.resize(numSamples);

    unsigned int sumKs = 0;
    std::vector <ScalarType> spherePoint(3,0.0);
    for (unsigned int i = 0;i < roundedNumLatitudes;++i)
    {
        double latitude = (i + 0.5) * M_PI / (2.0 * roundedNumLatitudes);
        unsigned int k5 = 0;

        if (i == roundedNumLatitudes - 1)
            k5 = numSamples - sumKs;
        else
            k5 = std::round(2.0 * numSamples * std::sin(latitude) * std::sin(M_PI / (4.0 * roundedNumLatitudes)));

        spherePoint[2] = std::cos(latitude);
        for (unsigned int j = 0;j < k5;++j)
        {
            double longitude = (j + 0.5) * 2.0 * M_PI / k5;

            spherePoint[0] = std::sin(latitude) * std::cos(longitude);
            spherePoint[1] = std::sin(latitude) * std::sin(longitude);

            spherePoints[sumKs + j] = spherePoint;
        }

        sumKs += k5;
    }
}

} //end namespace anima
