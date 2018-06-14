#pragma once

#include <vnl/vnl_vector_fixed.h>
#include <itkPoint.h>
#include <itkVector.h>

namespace anima
{
    
    template <class VectorType, class ScalarType>
    double EvaluateWatsonPDF(const VectorType &v, const VectorType &meanAxis, const ScalarType &kappa);
    
    template <class ScalarType>
    double EvaluateWatsonPDF(const vnl_vector_fixed <ScalarType,3> &v, const vnl_vector_fixed <ScalarType,3> &meanAxis, const ScalarType &kappa);
    
    template <class ScalarType>
    double EvaluateWatsonPDF(const itk::Point <ScalarType,3> &v, const itk::Point <ScalarType,3> &meanAxis, const ScalarType &kappa);
    
    template <class ScalarType>
    double EvaluateWatsonPDF(const itk::Vector <ScalarType,3> &v, const itk::Vector <ScalarType,3> &meanAxis, const ScalarType &kappa);
    
    template <class ScalarType>
    void GetStandardWatsonSHCoefficients(const ScalarType k, std::vector<ScalarType> &coefficients, std::vector<ScalarType> &derivatives);
    
} // end of namespace anima

#include "animaWatsonDistribution.hxx"
