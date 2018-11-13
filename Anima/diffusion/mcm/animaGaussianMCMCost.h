#pragma once

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>

#include <animaBaseMCMCost.h>

#include <animaMultiCompartmentModel.h>
#include <AnimaMCMExport.h>

namespace anima
{    

/**
 * @brief Class for computing marginal and profile costs and derivatives.
 * This is not thread safe at all so be sure to intantiate one per thread
 */
class ANIMAMCM_EXPORT GaussianMCMCost : public anima::BaseMCMCost
{
public:
    /** Standard class typedefs. */
    typedef GaussianMCMCost Self;
    typedef anima::BaseMCMCost Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(GaussianMCMCost, anima::BaseMCMCost)

    //! Get residual values for a given set of parameters, returns a vector of residuals
    MeasureType GetValues(const ParametersType &parameters) ITK_OVERRIDE;

    //! For the current set of parameters, compute the cost function value, requires GetValues to be called first
    double GetCurrentCostValue() ITK_OVERRIDE;

    //! Get residual derivatives for a given set of parameters, returns a matrix of residuals derivatives
    void GetDerivativeMatrix(const ParametersType &parameters, DerivativeMatrixType &derivative) ITK_OVERRIDE;

    //! Get cost function derivatives from a derivative matrix obtained from GetDerivativeMatrix
    void GetCurrentDerivative(DerivativeMatrixType &derivativeMatrix, DerivativeType &derivative) ITK_OVERRIDE;

    void SetMarginalEstimation(bool val) {m_MarginalEstimation = val;}

protected:
    GaussianMCMCost()
    {
        m_MarginalEstimation = false;
    }

    virtual ~GaussianMCMCost() {}

private:
    GaussianMCMCost(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    //! Switch to compute marginal or profile estimation
    bool m_MarginalEstimation;

    MeasureType m_Residuals;
};

} // end namespace anima
