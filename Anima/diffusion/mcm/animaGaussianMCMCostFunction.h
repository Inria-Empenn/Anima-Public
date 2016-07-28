#pragma once

#include <animaBaseMCMCostFunction.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>
#include <AnimaMCMExport.h>

namespace anima
{    

class ANIMAMCM_EXPORT GaussianMCMCostFunction : public BaseMCMCostFunction
{
public:
    /** Standard class typedefs. */
    typedef GaussianMCMCostFunction Self;
    typedef BaseMCMCostFunction Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(GaussianMCMCostFunction, BaseMCMCostFunction)

    typedef Superclass::MeasureType    MeasureType;
    typedef Superclass::DerivativeType DerivativeType;
    typedef Superclass::ParametersType ParametersType;

    // Object parameters contains the MCM parameters
    virtual MeasureType GetValue(const ParametersType &parameters) const ITK_OVERRIDE;
    virtual void GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const ITK_OVERRIDE;

    itkSetMacro(MarginalEstimation, bool)

protected:
    GaussianMCMCostFunction()
    {
        m_MarginalEstimation = true;
    }

    virtual ~GaussianMCMCostFunction() {}

    void GetMarginalDerivative(DerivativeType &derivative) const;
    void GetProfileDerivative(DerivativeType &derivative) const;

private:
    GaussianMCMCostFunction(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    //! Switch to compute marginal or profile estimation
    bool m_MarginalEstimation;

    // Utility variables to make ML estimation faster
    mutable double m_PredictedSquaredNorm;
};

} // end namespace anima
