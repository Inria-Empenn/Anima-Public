#pragma once
#include <cmath>

#include "animaWatsonDistribution.h"
#include <animaVectorOperations.h>
#include <animaErrorFunctions.h>

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
    
//    double b0 = boost::math::cyl_bessel_i(0, k / 2.0);
//    double b1 = boost::math::cyl_bessel_i(1, k / 2.0);
//    coefficients[0] = b0;
//    coefficients[1] = std::sqrt(5.0) / 4.0 * (b0 + 3.0 * b1);
//    coefficients[2] = 3.0 / 16.0 * (11.0 * b0 + 5.0 * (k - 7.0) * b1 / k);
//    coefficients[3] = std::sqrt(13.0) / (64.0 * k * k) * (11.0 * k * (2.0 * k - 21.0) * b0 + 21.0 * (44.0 + k * (2.0 * k - 3.0)) * b1);
//    coefficients[4] = std::sqrt(17.0) / (512.0 * k * k * k) * (k * (19305.0 - 858.0 * k + 326.0 * k * k) * b0 + 6.0 * (-12870.0 + k * (572.0 + k * (-594.0 + 31.0 * k))) * b1);
//    coefficients[5] = std::sqrt(21.0) / (2048.0 * k * k * k * k) * (k * (-1108536.0 + k * (36465.0 + 2.0 * k * (-9867.0 + 386.0 * k))) * b0 + 11.0 * (403104.0 + k * (-13260.0 + k * (19773.0 + 2.0 * k * (-325.0 + 58.0 * k)))) * b1);
//    coefficients[6] = 5.0 / (8192.0 * k * k * k * k * k) * (k * (81124680.0 + k * (-2116296.0 + k * (1625013.0 - 40664.0 * k + 5020.0 * k * k))) * b0 + 13.0 * (-24961440.0 + k * (651168.0 + k * (-1280049.0 + k * (32861.0 + 4.0 * k * (-2229.0 + 61.0 * k))))) * b1);
    
    double dawsonValue = anima::EvaluateDawsonIntegral(std::sqrt(k), true);
    double k2 = k * k;
    double k3 = k2 * k;
    double k4 = k3 * k;
    double k5 = k4 * k;
    double k6 = k5 * k;
    coefficients[0] = M_PI * dawsonValue;
    coefficients[1] = std::sqrt(5.0) * (6.0 * std::sqrt(k * M_PI) - (3.0 + 2.0 * k) * M_PI * dawsonValue) / 4.0;
    coefficients[2] = (30.0 * std::sqrt(k * M_PI) * (-21.0 + 2.0 * k) + 9.0 * (35.0 + 20.0 * k + 4.0 * k2) * M_PI * dawsonValue) / 32.0;
    coefficients[3] = std::sqrt(13.0) * (42.0 * std::sqrt(k * M_PI) * (165.0 - 20.0 * k + 4.0 * k2) - 5.0 * (693.0 + 378.0 * k + 84.0 * k2 + 8.0 * k3) * M_PI * dawsonValue) / 128.0;
    coefficients[4] = std::sqrt(17.0) * (6.0 * std::sqrt(k * M_PI) * (-225225.0 + 30030.0 * k - 7700.0 * k2 + 248.0 * k3) + 35.0 * (19305.0 + 10296.0 * k + 2376.0 * k2 + 288.0 * k3 + 16.0 * k4) * M_PI * dawsonValue) / 2048.0;
    coefficients[5] = std::sqrt(21.0) * (22.0 * std::sqrt(k * M_PI) * (3968055.0 - 556920.0 * k + 157248.0 * k2 - 7488.0 * k3 + 464.0 * k4) - 63.0 * (692835.0 + 364650.0 * k + 85800.0 * k2 + 11440.0 * k3 + 880.0 * k4 + 32.0 * k5) * M_PI * dawsonValue) / 8192.0;
    coefficients[6] = (130.0 * std::sqrt(k * M_PI) * (-540571185.0 + 78343650.0 * k - 23279256.0 * k2 + 1319472.0 * k3 - 119504.0 * k4 + 1952.0 * k5) + 1155.0 * (30421755.0 + 15872220.0 * k + 3779100.0 * k2 + 530400.0 * k3 + 46800.0 * k4 + 2496.0 * k5 + 64.0 * k6) * M_PI * dawsonValue) / 65536.0;
    
    for (unsigned int i = 0;i < nbCoefs;++i)
        coefficients[i] *= (std::exp(k) / std::pow(k, (2.0 * k + 1.0) / 2.0));
}

} // end namespace anima
