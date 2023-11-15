#pragma once

#include "AnimaSpecialFunctionsExport.h"

#include <cmath>

namespace anima
{

    ANIMASPECIALFUNCTIONS_EXPORT
    double PochHammer(const double &x,
                      const unsigned int n);

    //! According to Muller, K. E. (2001) ‘Computing the confluent hypergeometric function, M (a, b, x)’, Numerische Mathematik, pp. 179–196. Method 1.C, p.5
    ANIMASPECIALFUNCTIONS_EXPORT
    double
    KummerMethod1(const double &x, const double &a, const double &b,
                  const unsigned int maxIter = 1000, const double tol = 1.0e-8);

    //! According to Muller, K. E. (2001) ‘Computing the confluent hypergeometric function, M (a, b, x)’, Numerische Mathematik, pp. 179–196. Method 2, p.6
    ANIMASPECIALFUNCTIONS_EXPORT
    double
    KummerMethod2(const double &x, const double &a, const double &b,
                  const unsigned int maxIter = 1000, const double tol = 1.0e-8);

    //! According to Muller, K. E. (2001) ‘Computing the confluent hypergeometric function, M (a, b, x)’, Numerische Mathematik, pp. 179–196. Method with integral if b > a > 0
    ANIMASPECIALFUNCTIONS_EXPORT
    double
    KummerIntegrandMethod(const double &x, const double &a, const double &b);

    //! Computes the confluent hypergeometric function 1F1 also known as the Kummer function M. It calls Kummer function core to get the value.
    ANIMASPECIALFUNCTIONS_EXPORT
    double
    GetKummerFunctionValue(const double &x, const double &a, const double &b,
                           const unsigned int maxIter = 1000, const double tol = 1.0e-8);

    //! Computes the confluent hypergeometric function 1F1 also known as the Kummer function M. It returns a scaled value: exp(-x) * M(x,a,b).
    ANIMASPECIALFUNCTIONS_EXPORT
    double
    GetScaledKummerFunctionValue(const double &x, const double &a, const double &b,
                                 const unsigned int maxIter = 1000, const double tol = 1.0e-8);

    class KummerIntegrand
    {
    public:
        void SetXValue(double val) { m_XValue = val; }
        void SetAValue(double val) { m_AValue = val; }
        void SetBValue(double val) { m_BValue = val; }

        double operator()(const double t)
        {
            // For better numerical stability
            double tModified = t - (m_XValue > 0.0);
            return std::exp(m_XValue * tModified) * std::pow(t, m_AValue - 1.0) * std::pow(1.0 - t, m_BValue - m_AValue - 1.0);
        }

    private:
        double m_XValue, m_AValue, m_BValue;
    };

} // end namespace anima
