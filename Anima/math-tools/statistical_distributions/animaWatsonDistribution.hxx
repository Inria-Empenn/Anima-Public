#pragma once
#include <cmath>

#include "animaWatsonDistribution.h"
#include <animaVectorOperations.h>
#include <animaErrorFunctions.h>
#include <animaKummerFunctions.h>

#include <itkExceptionObject.h>
#include <itkObjectFactory.h>

#include <boost/math/special_functions/bessel.hpp>

namespace anima
{

template <class VectorType, class ScalarType>
double EvaluateWatsonPDF(const VectorType &v, const VectorType &meanAxis, const ScalarType &kappa)
{
    /************************************************************************************************
     * \fn template <class VectorType, class ScalarType>
     * 	   double
     *     EvaluateWatsonPDF(const VectorType &v,
     * 						 const VectorType &meanAxis,
     * 						 const ScalarType &kappa)
     *
     * \brief	Evaluate the Watson probability density function using the definition of
     * 			Fisher et al., Statistical Analysis of Spherical Data, 1993, p.89.
     *
     * \author	Aymeric Stamm
     * \date	July 2014
     *
     * \param	v				Sample on which evaluating the PDF.
     * \param	meanAxis        Mean axis of the Watson distribution.
     * \param	kappa		   	Concentration parameter of the Watson distribution.
     **************************************************************************************************/

    if (std::abs(kappa) < 1.0e-6)
        return 1.0 / (4.0 * M_PI);
    else if (kappa > 0)
    {
        double kappaSqrt = std::sqrt(kappa);
        double c = anima::ComputeScalarProduct(v, meanAxis);
        double inExp = kappa * (c * c - 1.0);
        return kappaSqrt * std::exp(inExp) / (4.0 * M_PI * anima::EvaluateDawsonFunctionNR(kappaSqrt));
    }
    else
    {
        double Ck = std::sqrt(-kappa / M_PI) / (2.0 * M_PI * std::erf(std::sqrt(-kappa)));
        double c = anima::ComputeScalarProduct(v, meanAxis);
        double inExp = kappa * c * c;
        return Ck * std::exp(inExp);
    }
}

template <class ScalarType>
double EvaluateWatsonPDF(const vnl_vector_fixed <ScalarType,3> &v, const vnl_vector_fixed <ScalarType,3> &meanAxis, const ScalarType &kappa)
{
    if (std::abs(v.squared_magnitude() - 1.0) > 1.0e-6 || std::abs(meanAxis.squared_magnitude() - 1.0) > 1.0e-6)
        throw itk::ExceptionObject(__FILE__, __LINE__,"The Watson distribution is on the 2-sphere.",ITK_LOCATION);

    return EvaluateWatsonPDF<vnl_vector_fixed <ScalarType,3>, ScalarType>(v,meanAxis,kappa);
}

template <class ScalarType>
double EvaluateWatsonPDF(const itk::Point <ScalarType,3> &v, const itk::Point <ScalarType,3> &meanAxis, const ScalarType &kappa)
{
    if (std::abs(v.GetVnlVector().squared_magnitude() - 1.0) > 1.0e-6 || std::abs(meanAxis.GetVnlVector().squared_magnitude() - 1.0) > 1.0e-6)
        throw itk::ExceptionObject(__FILE__, __LINE__,"The Watson distribution is on the 2-sphere.",ITK_LOCATION);

    return EvaluateWatsonPDF<itk::Point <ScalarType,3>, ScalarType>(v,meanAxis,kappa);
}

template <class ScalarType>
double EvaluateWatsonPDF(const itk::Vector <ScalarType,3> &v, const itk::Vector <ScalarType,3> &meanAxis, const ScalarType &kappa)
{
    if (std::abs(v.GetSquaredNorm() - 1.0) > 1.0e-6 || std::abs(meanAxis.GetSquaredNorm() - 1.0) > 1.0e-6)
        throw itk::ExceptionObject(__FILE__, __LINE__,"The Watson distribution is on the 2-sphere.",ITK_LOCATION);

    return EvaluateWatsonPDF<itk::Vector <ScalarType,3>, ScalarType>(v,meanAxis,kappa);
}

template <class ScalarType>
void GetStandardWatsonSHCoefficients(const ScalarType k, std::vector<ScalarType> &coefficients)
{
    // Computes the first 7 non-zero SH coefficients of the standard Watson distribution.
    const unsigned int nbCoefs = 7;
    coefficients.resize(nbCoefs);
    
    double kummerValue = anima::KummerFunction(k, 0.5, 1.5, true);
    double dawsonValue = M_PI * anima::EvaluateDawsonIntegral(std::sqrt(k), true);
    double kpiSqrt = std::sqrt(k * M_PI);
    double k2 = k * k;
    double k3 = k2 * k;
    double k4 = k3 * k;
    double k5 = k4 * k;
    double k6 = k5 * k;
    coefficients[0] = dawsonValue;
    coefficients[1] = (6.0 * kpiSqrt - (3.0 + 2.0 * k) * dawsonValue) / 4.0;
    coefficients[2] = (10.0 * kpiSqrt * (-21.0 + 2.0 * k) + 3.0 * (35.0 + 20.0 * k + 4.0 * k2) * dawsonValue) / 32.0;
    coefficients[3] = (42.0 * kpiSqrt * (165.0 - 20.0 * k + 4.0 * k2) - 5.0 * (693.0 + 378.0 * k + 84.0 * k2 + 8.0 * k3) * dawsonValue) / 128.0;
    coefficients[4] = (6.0 * kpiSqrt * (-225225.0 + 30030.0 * k - 7700.0 * k2 + 248.0 * k3) + 35.0 * (19305.0 + 10296.0 * k + 2376.0 * k2 + 288.0 * k3 + 16.0 * k4) * dawsonValue) / 2048.0;
    coefficients[5] = (22.0 * kpiSqrt * (3968055.0 - 556920.0 * k + 157248.0 * k2 - 7488.0 * k3 + 464.0 * k4) - 63.0 * (692835.0 + 364650.0 * k + 85800.0 * k2 + 11440.0 * k3 + 880.0 * k4 + 32.0 * k5) * dawsonValue) / 8192.0;
    coefficients[6] = (26.0 * kpiSqrt * (-540571185.0 + 78343650.0 * k - 23279256.0 * k2 + 1319472.0 * k3 - 119504.0 * k4 + 1952.0 * k5) + 231.0 * (30421755.0 + 15872220.0 * k + 3779100.0 * k2 + 530400.0 * k3 + 46800.0 * k4 + 2496.0 * k5 + 64.0 * k6) * dawsonValue) / 65536.0;
    
    for (unsigned int i = 0;i < nbCoefs;++i)
        coefficients[i] *= (std::sqrt(4.0 * i + 1.0) / kummerValue / std::pow(k, (2.0 * k + 1.0) / 2.0));
}

} // end namespace anima
