#pragma once

#include "AnimaArbSpecialFunctionsExport.h"

namespace anima
{

//! Computes Kummer's M functions (i.e. confluent hypergeometric function 1F1) using the ARB library
ANIMAARBSPECIALFUNCTIONS_EXPORT double GetKummerM(const double x, const double a, const double b);

} // end namespace anima
