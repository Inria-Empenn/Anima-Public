#pragma once

namespace anima
{

//! Gyromagnetic ratio (in 10^6 rad/s/T), from Nelson, J., Nuclear Magnetic Resonance Spectroscopy, Prentice Hall, Londres, 2003
const double DiffusionGyromagneticRatio = 267.513;

//! Given gyromagneitc ratio in 10^6 rad/s/T, gradient strength in mT/m and deltas in ms, computes b-value
inline double GetBValueFromAcquisitionParameters(double smallDelta, double largeDelta, double gradientStrength)
{
    double alpha = DiffusionGyromagneticRatio * smallDelta * gradientStrength;
    return alpha * alpha * (largeDelta - smallDelta / 3.0) * 1.0e-9;
}

inline double GetGradientStrengthFromBValue(double bValue, double smallDelta, double largeDelta)
{
    double gradientStrength = std::sqrt(bValue * 1.0e9 / (largeDelta - smallDelta / 3.0)) / (DiffusionGyromagneticRatio * smallDelta);
    return gradientStrength;
}

//! Default small delta value (classical values)
const double DiffusionSmallDelta = 10;

//! Default large delta value (classical values)
const double DiffusionLargeDelta = 40;

}
