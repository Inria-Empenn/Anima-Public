#pragma once
#include "animaSphereOperations.h"
#include <cmath>

namespace anima
{

template <class ScalarType> void GetSphereEvenSampling(std::vector < std::vector <ScalarType> > &spherePoints,
                                                       unsigned int numSamples)
{
    spherePoints.resize(numSamples);

    double offset = 2.0 / numSamples;
    double increment = M_PI * (3.0 - std::sqrt(5.0));

    std::vector <ScalarType> spherePoint(3,0.0);
    for (unsigned int i = 0;i < numSamples;++i)
    {
        spherePoint[1] = (i * offset) - 1 + offset / 2.0;
        double r = std::sqrt(1.0 - spherePoint[1] * spherePoint[1]);

        double phi = (i + 1) * increment;
        if (i == numSamples - 1)
            phi = 0;

        spherePoint[0] = std::cos(phi) * r;
        spherePoint[2] = std::sin(phi) * r;

        spherePoints[i] = spherePoint;
    }
}

} //end namespace anima
