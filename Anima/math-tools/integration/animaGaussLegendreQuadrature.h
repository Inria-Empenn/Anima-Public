#pragma once
#include "AnimaIntegrationExport.h"
#include <vector>

namespace anima
{

/**
 * @brief Computes the Gauss Laguerre quadrature of a function defined from R^+ into R.
 * Recenters the function on the interest zone with an affine relation, then uses Gauss Legendre on the left out part of the
 * function and computes the main part with Gauss Laguerre.
 */
class ANIMAINTEGRATION_EXPORT GaussLegendreQuadrature
{
public:
    GaussLegendreQuadrature();
    virtual ~GaussLegendreQuadrature() {}

    //! Specifies region on which the main part of the function is to be seen. If not specified, [-1,1] is the region
    void SetInterestZone(double minVal, double maxVal);

    void SetNumberOfComponents(unsigned int num) {m_NumberOfComponents = num;}

    template <class FunctionType>
    std::vector <double> GetVectorIntegralValue(FunctionType integrand)
    {
        std::vector <double> resVal(m_NumberOfComponents, 0.0);
        std::vector <double> tmpVal;

        for (unsigned int i = 0;i < m_Abcissas.size();++i)
        {
            tmpVal = integrand(m_Abcissas[i] * m_Slope + m_Intercept);
            for (unsigned int j = 0;j < m_NumberOfComponents;++j)
                resVal[j] += m_Slope * m_ValueWeights[i] * tmpVal[j];
        }

        return resVal;
    }

    template <class FunctionType>
    double GetIntegralValue(FunctionType integrand)
    {
        double resVal = 0.0;

        for (unsigned int i = 0;i < m_Abcissas.size();++i)
            resVal += m_Slope * m_ValueWeights[i] * integrand(m_Abcissas[i] * m_Slope + m_Intercept);

        return resVal;
    }

private:
    double m_Slope, m_Intercept;
    unsigned int m_NumberOfComponents;

    //! Function value weights as computed in Lowan, A. et al. (1942).
    const std::vector <double> m_ValueWeights =
    {
        0.030753241996117,
        0.070366047488108,
        0.107159220467172,
        0.139570677926154,
        0.166269205816994,
        0.186161000015562,
        0.198431485327111,
        0.202578241925561,
        0.198431485327111,
        0.186161000015562,
        0.166269205816994,
        0.139570677926154,
        0.107159220467172,
        0.070366047488108,
        0.030753241996117
    };

    //! Zeroes of Legendre polynomials of order 15 as computed in Lowan, A. et al. (1942).
    const std::vector <double> m_Abcissas =
    {
        -0.987992518020485,
        -0.937273392400706,
        -0.848206583410427,
        -0.724417731360170,
        -0.570972172608539,
        -0.394151347077563,
        -0.201194093997435,
        0.0,
        0.201194093997435,
        0.394151347077563,
        0.570972172608539,
        0.724417731360170,
        0.848206583410427,
        0.937273392400706,
        0.987992518020485
    };
};

} // end of namespace anima
