/*
 *  animaKummerFunction.h
 *
 *
 *  Created by Aymeric Stamm on 07/12/2017.
 *  Copyright 2017 CRL. All rights reserved.
 *
 */

#pragma once

#include "AnimaSpecialFunctionsExport.h"

namespace anima
{

ANIMASPECIALFUNCTIONS_EXPORT
double PochHammer(const double &x,
                  const unsigned int n);

ANIMASPECIALFUNCTIONS_EXPORT
double
KummerMethod1(const double &x,
              const double &a,
              const double &b,
              const unsigned int maxIter = 10000,
              const double tol = 1.0e-15);

ANIMASPECIALFUNCTIONS_EXPORT
double
KummerMethod2(const double &x,
              const double &a,
              const double &b,
              const unsigned int maxIter = 10000,
              const double tol = 1.0e-15);

ANIMASPECIALFUNCTIONS_EXPORT
double
KummerFunction(const double &x,
               const double &a,
               const double &b,
               const unsigned int maxIter = 10000,
               const double tol = 1.0e-15);

} // end namespace anima
