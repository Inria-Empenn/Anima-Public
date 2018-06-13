#pragma once

namespace anima
{

//! Gyromagnetic ratio (in rad/ms/T), from Nelson, J., Nuclear Magnetic Resonance Spectroscopy, Prentice Hall, Londres, 2003
const double DiffusionGyromagneticRatio = 267513.0;

//! Given gyromagnetic ratio in rad/ms/T, gradient strength in T/mm and deltas in ms, computes b-value in s/mm^2
inline double GetBValueFromAcquisitionParameters(double smallDelta, double largeDelta, double gradientStrength)
{
    double alpha = DiffusionGyromagneticRatio * smallDelta * gradientStrength;
    return alpha * alpha * (largeDelta - smallDelta / 3.0) * 1.0e-3;
}

//! Given b-value in s/mm^2 and deltas in ms, computes gradient strength in T/mm
inline double GetGradientStrengthFromBValue(double bValue, double smallDelta, double largeDelta)
{
    double alpha = std::sqrt(bValue * 1.0e3 / (largeDelta - smallDelta / 3.0));
    return alpha / (DiffusionGyromagneticRatio * smallDelta);
}

//! Default small delta value (classical values)
const double DiffusionSmallDelta = 10.0;

//! Default large delta value (classical values)
const double DiffusionLargeDelta = 40.0;

}
