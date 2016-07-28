#pragma once

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>

#include <itkLightObject.h>
#include <itkArray2D.h>
#include <itkOptimizerParameters.h>

#include <animaMultiCompartmentModel.h>
#include <AnimaMCMExport.h>

namespace anima
{    

/**
 * @brief Class for computing variable projection costs and derivatives. Right now, it is only available for Gaussian noise.
 * In the future, it could be extended to handle other noise configurations by creating a base class from this one. By the way, this is not thread safe at all
 * so be sure to intantiate one per thread
 */
class ANIMAMCM_EXPORT GaussianMCMVariableProjectionCost : public itk::LightObject
{
public:
    /** Standard class typedefs. */
    typedef GaussianMCMVariableProjectionCost Self;
    typedef itk::LightObject Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(GaussianMCMVariableProjectionCost, itk::LightObject)

    typedef itk::Array2D <double> DerivativeMatrixType;
    typedef itk::Array <double> MeasureType;
    typedef itk::Array <double> DerivativeType;
    typedef itk::OptimizerParameters <double> ParametersType;

    typedef anima::MultiCompartmentModel MCMType;
    typedef MCMType::Pointer MCMPointer;
    typedef MCMType::Vector3DType Vector3DType;
    typedef MCMType::ListType ListType;

    void SetObservedSignals(ListType &value) {m_ObservedSignals = value;}
    void SetGradients(std::vector<Vector3DType> &value) {m_Gradients = value;}
    void SetBValues(ListType &value) {m_BValues = value;}

    void SetMCMStructure(MCMType *model) {m_MCMStructure = model;}
    MCMPointer &GetMCMStructure() {return m_MCMStructure;}

    //! Get residual values for a given set of parameters, returns a vector of residuals
    MeasureType GetValues(const ParametersType &parameters);

    //! For the current set of parameters, compute the cost function value, requires GetValues to be called first
    double GetCurrentCostValue();

    //! Get residual derivatives for a given set of parameters, returns a matrix of residuals derivatives
    void GetDerivativeMatrix(const ParametersType &parameters, DerivativeMatrixType &derivative);

    //! Get cost function derivatives from a derivative matrix obtained from GetDerivativeMatrix
    void GetCurrentDerivative(DerivativeMatrixType &derivativeMatrix, DerivativeType &derivative);

    //! Returns number of optimized parameters
    unsigned int GetNumberOfParameters() const
    {
        if (!m_MCMStructure)
            itkExceptionMacro("A multi-compartment model has to be set, no computation can be done otherwise.");

        return m_MCMStructure->GetNumberOfParameters();
    }

    unsigned int GetNumberOfObservations() const
    {
        return m_ObservedSignals.size();
    }

    double GetB0Value() {return m_B0Value;}
    std::vector <double> &GetOptimalWeights() {return m_OptimalWeights;}
    double GetSigmaSquare() {return m_SigmaSquare;}

    void SetUseDerivative(bool derivative) {m_UseDerivative = derivative;}

protected:
    GaussianMCMVariableProjectionCost()
    {
        m_B0Value = 0;
        m_SigmaSquare = 0;
        m_UseDerivative = false;
    }

    virtual ~GaussianMCMVariableProjectionCost() {}

    //! Computes maximum likelihood estimates of weights and B0
    void SolveUnconstrainedLinearLeastSquares();

    //! Computes maximum likelihood estimates of weights and B0 on the borders of the domain (has to be used after SolveMaximumLikelihoodPDE)
    void SolveLinearLeastSquaresOnBorders();

    bool CheckBoundaryConditions();
    void PrepareDataForLLS();

private:
    GaussianMCMVariableProjectionCost(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    double m_B0Value;
    double m_SigmaSquare;
    std::vector <double> m_PredictedSignals;

    ListType m_TestedParameters;

    ListType m_ObservedSignals;
    std::vector<Vector3DType> m_Gradients;
    ListType m_BValues;

    MCMPointer m_MCMStructure;
    bool m_UseDerivative;

    // Utility variables to make ML estimation faster
    vnl_matrix<double> m_PhiInverseG, m_PhiInverseGCopy;
    std::vector <double> m_OptimalWeights, m_OptimalWeightsCopy;
    MeasureType m_Residuals, m_ResidualsCopy;
    std::vector <bool> m_CompartmentSwitches, m_CompartmentSwitchesCopy;
    std::vector <double> m_ProjFirstCompartmentSignals, m_ProjFirstCompartmentSignalsCopy;
    
    std::vector <unsigned int> m_IndexesUsefulCompartments;
    std::vector < std::vector <double> > m_PredictedSignalAttenuations;
    std::vector< vnl_matrix<double> > m_SignalAttenuationsJacobian;
    double m_ProjFirstCompartmentSqNorm;
};

} // end namespace anima
