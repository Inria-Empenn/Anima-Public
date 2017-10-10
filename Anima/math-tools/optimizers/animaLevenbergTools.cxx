#include <animaLevenbergTools.h>
#include <cmath>

namespace anima
{

namespace levenberg
{

double UnboundValue(double x, double lowerBound, double upperBound)
{
    if (x < lowerBound)
        x = lowerBound;
    if (x > upperBound)
        x = upperBound;

    return std::asin(2.0 * (x - lowerBound) / (upperBound - lowerBound) - 1.0);
}

double ComputeBoundedValue(double x, double &inputSign, double lowerBound, double upperBound)
{
    double sinVal = std::sin(x);

    if (std::cos(x) >= 0)
        inputSign = 1.0;
    else
        inputSign = -1.0;

    return (upperBound - lowerBound) * (sinVal + 1.0) / 2.0 + lowerBound;
}

double BoundedDerivativeAddOn(double x, double inputSign, double lowerBound, double upperBound)
{
    if (x < lowerBound)
        x = lowerBound;
    if (x > upperBound)
        x = upperBound;

    return inputSign * std::sqrt((x - lowerBound) * (upperBound - x));
}

}

}
