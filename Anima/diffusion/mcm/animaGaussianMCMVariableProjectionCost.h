#pragma once

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>

#include <animaBaseMCMCost.h>
#include <animaNNLSOptimizer.h>
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

protected:
    GaussianMCMVariableProjectionCost()
    {
        m_NNLSBordersOptimizer = anima::NNLSOptimizer::New();
    }

    virtual ~GaussianMCMVariableProjectionCost() ITK_OVERRIDE {}

    //! Computes maximum likelihood estimates of weights
    void SolveLinearLeastSquares();

    bool CheckBoundaryConditions();
    void PrepareDataForLLS();
    void PrepareDataForDerivative();

private:
    GaussianMCMVariableProjectionCost(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    // Utility variables to make ML estimation faster
    vnl_matrix <double> m_FMatrixInverseG;
    std::vector <double> m_OptimalWeights;
    MeasureType m_Residuals;
    std::vector <bool> m_CompartmentSwitches;
    ListType m_FSignal;
    vnl_matrix <double> m_GramMatrix, m_InverseGramMatrix;
    ParametersType m_VNLSignals, m_OptimalNNLSWeights;
    anima::NNLSOptimizer::Pointer m_NNLSBordersOptimizer;
    
    std::vector <unsigned int> m_IndexesUsefulCompartments;
    vnl_matrix <double> m_PredictedSignalAttenuations;
    std::vector< vnl_matrix<double> > m_SignalAttenuationsJacobian;
};

} // end namespace anima
