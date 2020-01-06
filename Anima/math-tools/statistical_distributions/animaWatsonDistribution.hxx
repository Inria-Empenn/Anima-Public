#pragma once
#include <cmath>

#include "animaWatsonDistribution.h"
#include <animaVectorOperations.h>
#include <animaErrorFunctions.h>

#include <itkMacro.h>
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
    void GetStandardWatsonSHCoefficients(const ScalarType k, std::vector<ScalarType> &coefficients, std::vector<ScalarType> &derivatives)
{
    // Computes the first 7 non-zero SH coefficients of the standard Watson PDF (multiplied by 4 M_PI).
    const unsigned int nbCoefs = 7;
    coefficients.resize(nbCoefs);
    derivatives.resize(nbCoefs);
    
    double sqrtPi = std::sqrt(M_PI);
    double dawsonValue = anima::EvaluateDawsonIntegral(std::sqrt(k), true);
    double k2 = k * k;
    double k3 = k2 * k;
    double k4 = k3 * k;
    double k5 = k4 * k;
    double k6 = k5 * k;
    double k7 = k6 * k;
    
    coefficients[0] = 2.0 * sqrtPi;
    derivatives[0] = 0.0;
    
    if (k < 1.0e-4) // Handles small k values by Taylor series expansion
    {
        coefficients[1] = k / 3.0;
        coefficients[2] = k2 / 35.0;
        coefficients[3] = k3 / 693.0;
        coefficients[4] = k4 / 19305.0;
        coefficients[5] = k5 / 692835.0;
        coefficients[6] = k6 / 30421755.0;
        derivatives[1] = 1.0 / 3.0;
        derivatives[2] = 2.0 * k / 35.0;
        derivatives[3] = k2 / 231.0;
        derivatives[4] = 4.0 * k3 / 19305.0;
        derivatives[5] = k4 / 138567.0;
        derivatives[6] = 6.0 * k5 / 30421755.0;
        
        for (unsigned int i = 1;i < nbCoefs;++i)
        {
            double tmpVal = std::pow(2.0, i + 1.0) * sqrtPi / std::sqrt(4.0 * i + 1.0);
            coefficients[i] *= tmpVal;
            derivatives[i] *= tmpVal;
        }
        
        return;
    }
    
    coefficients[1] = (3.0 - (3.0 + 2.0 * k) * dawsonValue) / (2.0 * k);
    coefficients[2] = (5.0 * (-21.0 + 2.0 * k) + 3.0 * (35.0 + 4.0 * k * (5.0 + k)) * dawsonValue) / (16.0 * k2);
    coefficients[3] = (21.0 * (165.0 + 4.0 * (-5.0 + k) * k) - 5.0 * (693.0 + 378.0 * k + 84.0 * k2 + 8.0 * k3) * dawsonValue) / (64.0 * k3);
    coefficients[4] = (3.0 * (-225225.0 + 2.0 * k * (15015.0 + 2.0 * k * (-1925.0 + 62.0 * k))) + 35.0 * (19305.0 + 8.0 * k * (1287.0 + k * (297.0 + 2.0 * k *  (18.0 + k)))) * dawsonValue) / (1024.0 * k4);
    coefficients[5] = (11.0 * (3968055.0 + 8.0 * k * (-69615.0 + 2.0 * k * (9828.0 + k * (-468.0 + 29.0 * k)))) - 63.0 * (692835.0 + 2.0 * k * (182325.0 + 4.0 * k * (10725.0 + 2.0 * k * (715.0 + k * (55.0 + 2.0 * k))))) * dawsonValue) / (4096.0 * k5);
    coefficients[6] = (13.0 * (-540571185.0 + 2.0 * k * (39171825.0 + 4.0 * k * (-2909907.0 + 2.0 * k * (82467.0 + k * (-7469.0 + 122.0 * k))))) + 231.0 * (30421755.0 + 4.0 * k * (3968055.0 + k * (944775.0 + 4.0 * k * (33150.0 + k * (2925.0 + 4.0 * k * (39.0 + k)))))) * dawsonValue) / (32768.0 * k6);
    
    derivatives[1] = 3.0 * (-1.0 + (-1.0 + 2.0 * k) * dawsonValue + 2.0 * dawsonValue * dawsonValue) / (4.0 * k2);
    derivatives[2] = 5.0 * ((21.0 - 2.0 * k) + (63.0 + 4.0 * (-11.0 + k) * k) * dawsonValue - 12.0 * (7.0 + 2.0 * k) * dawsonValue * dawsonValue) / (32.0 * k3);
    derivatives[3] = 21.0 * ((-165.0 - 4.0 * (-5.0 + k) * k) + (-5.0 + 2.0 * k) * (165.0 + 4.0 * (-3.0 + k) * k) * dawsonValue + 10.0 * (99.0 + 4.0 * k * (9.0 + k)) * dawsonValue * dawsonValue) / (128.0 * k4);
    derivatives[4] = 3.0 * ((225225.0 - 2.0 * k * (15015.0 + 2.0 * k * (-1925.0 + 62.0 * k))) + (1576575.0 + 8.0 * k * (-75075.0 + k * (10395.0 + 2.0 * k * (-978.0 + 31.0 * k)))) * dawsonValue - 840.0 * (2145.0 + 2.0 * k * (429.0 + 66.0 * k + 4.0 * k2)) * dawsonValue * dawsonValue) / (2048.0 * k5);
    derivatives[5] = 11.0 * ((-3968055.0 - 8.0 * k * (-69615.0 + 2.0 * k * (9828.0 + k * (-468.0 + 29.0 * k)))) + (-35712495.0 + 2.0 * k * (5917275.0 + 8.0 * k * (-118755.0 + k * (21060.0 + k * (-965.0 + 58.0 * k))))) * dawsonValue + 630.0 * (62985.0 + 8.0 * k * (3315.0 + k * (585.0 + 2.0 * k * (26.0 + k)))) * dawsonValue * dawsonValue) / (8192.0 * k6);
    derivatives[6] = 13.0 * ((540571185.0 - 2.0 * k * (39171825.0 + 4.0 * k * (-2909907.0 + 2.0 * k * (82467.0 + k * (-7469.0 + 122.0 * k))))) + (5946283035.0 + 4.0 * k * (-446558805.0 + k * (79910523.0 + 4.0 * k * (-3322242.0 + k * (187341.0 + 4.0 * k * (-3765.0 + 61.0 * k)))))) * dawsonValue - 2772.0 * (2340135.0 + 2.0 * k * (508725.0 + 4.0 * k * (24225.0 + 2.0 * k * (1275.0 + k * (75.0 + 2.0 * k))))) * dawsonValue * dawsonValue) / (65536.0 * k7);
    
    for (unsigned int i = 1;i < nbCoefs;++i)
    {
        double sqrtVal = std::sqrt(1.0 + 4.0 * i);
        coefficients[i] *= (sqrtPi * sqrtVal / dawsonValue);
        derivatives[i] *= (sqrtPi * sqrtVal / (dawsonValue * dawsonValue));
    }
}

} // end namespace anima
