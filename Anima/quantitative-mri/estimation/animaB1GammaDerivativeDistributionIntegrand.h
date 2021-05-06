#pragma once

#include "AnimaRelaxometryExport.h"

#include <vector>
#include <map>

#include <animaB1GammaDistributionIntegrand.h>

namespace anima
{

/**
 * \class B1GammaDerivativeDistributionIntegrand
 * @brief Integrand to compute the internal derivative integral per distribution in B1GammaMixtureT2RelaxometryCostFunction
 */
class ANIMARELAXOMETRY_EXPORT B1GammaDerivativeDistributionIntegrand : public B1GammaDistributionIntegrand
{
public:
    using Superclass = B1GammaDistributionIntegrand;
    B1GammaDerivativeDistributionIntegrand()
    {
        m_B1DerivativeFlag = true;
    }

    void SetB1DerivativeFlag(bool val) {m_B1DerivativeFlag = val;}

    virtual std::vector <double> operator() (double const t) override;

private:
    bool m_B1DerivativeFlag;
};

} // end namespace anima
