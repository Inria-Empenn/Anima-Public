#pragma once
#include <cmath>

#include "animaWatsonDistribution.h"
#include <animaVectorOperations.h>
#include <animaErrorFunctions.h>

#include <itkExceptionObject.h>
#include <itkObjectFactory.h>

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

} // end namespace anima
