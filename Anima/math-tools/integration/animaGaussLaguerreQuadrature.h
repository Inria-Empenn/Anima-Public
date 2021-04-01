#pragma once
#include "AnimaIntegrationExport.h"
#include <vector>

#include <animaGaussLegendreQuadrature.h>

namespace anima
{

/**
 * @brief Computes the Gauss Laguerre quadrature of a function defined from R^+ into R.
 * Recenters the function on the interest zone with an affine relation, then uses Gauss Legendre on the left out part of the
 * function and computes the main part with Gauss Laguerre.
 */
class ANIMAINTEGRATION_EXPORT GaussLaguerreQuadrature
{
public:
    GaussLaguerreQuadrature();
    virtual ~GaussLaguerreQuadrature() {}

    //! Specifies region on which the main part of the function is to be seen. If not specified, R^+ is the region
    void SetInterestZone(double minVal, double maxVal);

    void SetNumberOfComponents(unsigned int num) {m_NumberOfComponents = num;}

    template <class FunctionType>
    double GetIntegralValue(FunctionType integrand)
    {
        double deltaPart = 0.0;
        if (m_DeltaValue > 0.0)
        {
            anima::GaussLegendreQuadrature glQuad;
            glQuad.SetInterestZone(0.0, m_DeltaValue);
            deltaPart = glQuad.GetIntegralValue(integrand);
        }

        double laguerrePart = 0.0;
        for (unsigned int i = 0;i < m_Abcissas.size();++i)
            laguerrePart += m_ValueWeights[i] * integrand(m_ScaleValue * m_Abcissas[i] + m_DeltaValue) * m_ScaleValue;

        return deltaPart + laguerrePart;

    }

    template <class FunctionType>
    std::vector <double> GetVectorIntegralValue(FunctionType integrand)
    {
        std::vector <double> integrandPart;
        std::vector <double> resVal(m_NumberOfComponents, 0.0);

        if (m_DeltaValue > 0.0)
        {
            anima::GaussLegendreQuadrature glQuad;
            glQuad.SetInterestZone(0.0, m_DeltaValue);
            glQuad.SetNumberOfComponents(m_NumberOfComponents);
            resVal = glQuad.GetVectorIntegralValue(integrand);
        }

        for (unsigned int i = 0;i < m_Abcissas.size();++i)
        {
            integrandPart = integrand(m_ScaleValue * m_Abcissas[i] + m_DeltaValue);
            for (unsigned int j = 0;j < m_NumberOfComponents;++j)
                resVal[j] += m_ValueWeights[i] * integrandPart[j] * m_ScaleValue;
        }

        return resVal;
    }

private:
    //! Delta for computing the integral. Since the Laguerre zeroes are up to about 50, this is to ensure the function values are meaningful
    double m_DeltaValue;
    double m_ScaleValue;

    unsigned int m_NumberOfComponents;

    //! Function value weights as computed in Rabinowitz, P.; Weiss, G. (1959). doi:10.1090/S0025-5718-1959-0107992-3
    const std::vector <double> m_ValueWeights =
    {
        2.25030314864247252e-1,
        5.25836052762342454e-1,
        8.31961391687087088e-1,
        1.14609924096375170,
        1.47175131696680859,
        1.81313408738134816,
        2.17551751969460745,
        2.56576275016502921,
        2.99321508637137516,
        3.47123448310209029,
        4.02004408644466886,
        4.67251660773285420,
        5.48742065798615247,
        6.58536123328921366,
        8.27635798436423448,
        1.18242775516584348e1
    };

    //! Zeroes of Laguerre polynomials as computed in Rabinowitz, P.; Weiss, G. (1959). doi:10.1090/S0025-5718-1959-0107992-3
    const std::vector <double> m_Abcissas =
    {
        8.76494104789278403e-2,
        4.62696328915080832e-1,
        1.14105777483122686,
        2.12928364509838062,
        3.43708663389320665,
        5.07801861454976791,
        7.07033853504823413,
        9.43831433639193878,
        1.22142233688661587e1,
        1.54415273687816171e1,
        1.91801568567531349e1,
        2.35159056939919085e1,
        2.85787297428821404e1,
        3.45833987022866258e1,
        4.19404526476883326e1,
        5.17011603395433184e1
    };
};

} // end of namespace anima
