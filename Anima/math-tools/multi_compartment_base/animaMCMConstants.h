#pragma once

namespace anima
{

// The following units are in place to keep things right: distances are in millimeters, time in seconds, magnetic field in Tesla

//! Gyromagnetic ratio (in rad/s/T), from Nelson, J., Nuclear Magnetic Resonance Spectroscopy, Prentice Hall, Londres, 2003
const double DiffusionGyromagneticRatio = 267.513e6;

//! Given gyromagnetic ratio in rad/s/T, gradient strength in T/mm and deltas in s, computes b-value in s/mm^2
inline double GetBValueFromAcquisitionParameters(double smallDelta, double bigDelta, double gradientStrength)
{
    double alpha = DiffusionGyromagneticRatio * smallDelta * gradientStrength;
    return alpha * alpha * (bigDelta - smallDelta / 3.0);
}

//! Given b-value in s/mm^2 and deltas in s, computes gradient strength in T/mm
inline double GetGradientStrengthFromBValue(double bValue, double smallDelta, double bigDelta)
{
    double alpha = std::sqrt(bValue / (bigDelta - smallDelta / 3.0));
    return alpha / (DiffusionGyromagneticRatio * smallDelta);
}

//! Default small delta value (classical values)
const double DiffusionSmallDelta = 10.0e-3;

//! Default big delta value (classical values)
const double DiffusionBigDelta = 40.0e-3;

}
