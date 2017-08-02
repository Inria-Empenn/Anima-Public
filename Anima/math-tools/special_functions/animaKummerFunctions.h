#pragma once

#include "AnimaSpecialFunctionsExport.h"

namespace anima
{

ANIMASPECIALFUNCTIONS_EXPORT
double PochHammer(const double &x,
                  const unsigned int n);

//! According to Muller, K. E. (2001) ‘Computing the confluent hypergeometric function, M (a, b, x)’, Numerische Mathematik, pp. 179–196. Method 1, p.5
ANIMASPECIALFUNCTIONS_EXPORT
double
KummerMethod1(const double &x,
              const double &a,
              const double &b,
              const unsigned int maxIter = 10000,
              const double tol = 1.0e-15);

//! According to Muller, K. E. (2001) ‘Computing the confluent hypergeometric function, M (a, b, x)’, Numerische Mathematik, pp. 179–196. Method 2, p.6
ANIMASPECIALFUNCTIONS_EXPORT
double
KummerMethod2(const double &x,
              const double &a,
              const double &b,
              const unsigned int maxIter = 10000,
              const double tol = 1.0e-15);

//! According to Muller, K. E. (2001) ‘Computing the confluent hypergeometric function, M (a, b, x)’, Numerische Mathematik, pp. 179–196. Switch between Method 1 and 2 according to recommandation at the end of page 5. Covers most situations (at least all common situations encountered in MR image processing).
ANIMASPECIALFUNCTIONS_EXPORT
double
KummerFunction(const double &x,
               const double &a,
               const double &b,
               const unsigned int maxIter = 10000,
               const double tol = 1.0e-15);

} // end namespace anima
