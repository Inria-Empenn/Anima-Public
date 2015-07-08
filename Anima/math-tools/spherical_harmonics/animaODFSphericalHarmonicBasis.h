#pragma once

#include <animaSphericalHarmonic.h>
#include <itkVariableLengthVector.h>
#include <vector>
#include <AnimaSHToolsExport.h>

namespace anima
{
class ANIMASHTOOLS_EXPORT ODFSphericalHarmonicBasis
{
public:
    ODFSphericalHarmonicBasis(unsigned int L);
    virtual ~ODFSphericalHarmonicBasis() {m_SphericalHarmonics.clear();}

    // T has to be a vector type with the [] operator
    template <class T> double getValueAtPosition(const T &coefficients, double theta, double phi);

    template <class T> double getThetaFirstDerivativeValueAtPosition(const T &coefficients,
                                                  double theta, double phi);

    template <class T> double getPhiFirstDerivativeValueAtPosition(const T &coefficients,
                                                double theta, double phi);

    template <class T> double getThetaSecondDerivativeValueAtPosition(const T &coefficients,
                                                   double theta, double phi);

    template <class T> double getThetaPhiDerivativeValueAtPosition(const T &coefficients,
                                                double theta, double phi);

    template <class T> double getPhiSecondDerivativeValueAtPosition(const T &coefficients,
                                                 double theta, double phi);

    template <class T> double getCurvatureAtPosition(const T &coefficients,
                                  double theta, double phi);

    double getNthSHValueAtPosition(int k, int m, double theta, double phi);

    template <class T> itk::VariableLengthVector <T>
    GetSampleValues(itk::VariableLengthVector <T> &data,
                    std::vector < std::vector <double> > &m_SampleDirections);
private:
    unsigned int m_LOrder;
    std::vector < SphericalHarmonic > m_SphericalHarmonics;
};

} // end namespace odf

#include "animaODFSphericalHarmonicBasis.hxx"
