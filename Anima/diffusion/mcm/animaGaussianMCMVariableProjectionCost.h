#pragma once

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>

#include <animaBaseMCMCost.h>
#include <AnimaMCMExport.h>

namespace anima
{    

/**
 * @brief Class for computing variable projection costs and derivatives. Right now, it is only available for Gaussian noise.
 * By the way, this is not thread safe at all so be sure to intantiate one per thread
 */
class ANIMAMCM_EXPORT GaussianMCMVariableProjectionCost : public anima::BaseMCMCost
{
public:
    /** Standard class typedefs. */
    typedef GaussianMCMVariableProjectionCost Self;
    typedef anima::BaseMCMCost Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(GaussianMCMVariableProjectionCost, anima::BaseMCMCost)

    //! Get residual values for a given set of parameters, returns a vector of residuals
    MeasureType GetValues(const ParametersType &parameters) ITK_OVERRIDE;

    //! For the current set of parameters, compute the cost function value, requires GetValues to be called first
    double GetCurrentCostValue() ITK_OVERRIDE;

    //! Get residual derivatives for a given set of parameters, returns a matrix of residuals derivatives
    void GetDerivativeMatrix(const ParametersType &parameters, DerivativeMatrixType &derivative) ITK_OVERRIDE;

    //! Get cost function derivatives from a derivative matrix obtained from GetDerivativeMatrix
    void GetCurrentDerivative(DerivativeMatrixType &derivativeMatrix, DerivativeType &derivative) ITK_OVERRIDE;

    std::vector <double> &GetOptimalWeights() {return m_OptimalWeights;}
    void SetUseDerivative(bool derivative) {m_UseDerivative = derivative;}

protected:
    GaussianMCMVariableProjectionCost()
    {
        m_UseDerivative = false;
    }

    virtual ~GaussianMCMVariableProjectionCost() ITK_OVERRIDE {}

    //! Computes maximum likelihood estimates of weights
    void SolveUnconstrainedLinearLeastSquares();

    //! Computes maximum likelihood estimates of weights on the borders of the domain (has to be used after SolveUnconstrainedLinearLeastSquares)
    void SolveLinearLeastSquaresOnBorders();

    bool CheckBoundaryConditions();
    void PrepareDataForLLS();

private:
    GaussianMCMVariableProjectionCost(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    bool m_UseDerivative;

    // Utility variables to make ML estimation faster
    vnl_matrix <double> m_FMatrixInverseG, m_FMatrixInverseGCopy;
    std::vector <double> m_OptimalWeights, m_OptimalWeightsCopy;
    MeasureType m_Residuals, m_ResidualsCopy;
    std::vector <bool> m_CompartmentSwitches, m_CompartmentSwitchesCopy;
    ListType m_FSignal;
    vnl_matrix <double> m_CompleteGramMatrix, m_GramMatrix, m_InverseGramMatrix;
    
    std::vector <unsigned int> m_IndexesUsefulCompartments;
    std::vector < std::vector <double> > m_PredictedSignalAttenuations;
    std::vector< vnl_matrix<double> > m_SignalAttenuationsJacobian;
};

} // end namespace anima
