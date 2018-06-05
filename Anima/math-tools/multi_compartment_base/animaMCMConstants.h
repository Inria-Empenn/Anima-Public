#pragma once

namespace anima
{

//! Gyromgnetic ratio (in 10^6 rad/s/T), from Nelson, J., Nuclear Magnetic Resonance Spectroscopy, Prentice Hall, Londres, 2003
const double DiffusionGyromagneticRatio = 267.513;

inline double GetBValueFromAcquisitionParameters(double smallDelta, double largeDelta, double gradientStrength)
{
    double alpha = DiffusionGyromagneticRatio * smallDelta * gradientStrength;
    return alpha * alpha * (largeDelta - smallDelta / 3.0);
}

inline double GetGradientStrengthFromBValue(double bValue, double smallDelta, double largeDelta)
{
    double gradientStrength = std::sqrt(largeDelta - smallDelta / 3.0) / (smallDelta * DiffusionGyromagneticRatio);
    return gradientStrength;
}

//! Default small delta value to get bval = G^2
const double DiffusionSmallDelta = 1.0 / DiffusionGyromagneticRatio;

//! Default large delta value to get bval = G^2
const double DiffusionLargeDelta = 1.0 + DiffusionSmallDelta / 3.0;

}
