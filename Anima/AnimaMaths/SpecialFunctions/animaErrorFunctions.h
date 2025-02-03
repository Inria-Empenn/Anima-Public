#pragma once

#include <libAnimaMathsExport.h>

#include <cmath>

namespace anima
{

    class DawsonIntegrand
    {
    public:
        void SetXValue(double val) { m_XValue = val; }

        double operator()(const double t)
        {
            double inexpValue = m_XValue * m_XValue * (t * t - 1.0);
            return std::exp(inexpValue);
        }

    private:
        double m_XValue;
    };

    LIBANIMAMATHS_EXPORT double EvaluateDawsonIntegral(const double x, const bool scaled = false);

    LIBANIMAMATHS_EXPORT double EvaluateDawsonFunctionNR(double x);

    LIBANIMAMATHS_EXPORT double EvaluateDawsonFunction(double x);

    LIBANIMAMATHS_EXPORT double EvaluateWImFunction(double x);

    LIBANIMAMATHS_EXPORT double EvaluateWImY100Function(double y100, double x);

} // end namespace anima
