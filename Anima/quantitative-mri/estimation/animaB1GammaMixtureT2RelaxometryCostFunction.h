#pragma once

#include "AnimaRelaxometryExport.h"

#include <vector>
#include <itkSingleValuedCostFunction.h>
#include <vnl/vnl_matrix.h>

#include <animaEPGSignalSimulator.h>
#include <animaCholeskyDecomposition.h>
#include <animaNNLSOptimizer.h>

namespace anima
{

class ANIMARELAXOMETRY_EXPORT B1GammaMixtureT2RelaxometryCostFunction :
        public itk::SingleValuedCostFunction
{
public:
    /** Standard class typedefs. */
    typedef B1GammaMixtureT2RelaxometryCostFunction Self;
    typedef SingleValuedCostFunction Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(B1GammaMixtureT2RelaxometryCostFunction, SingleValuedCostFunction)

    typedef Superclass::MeasureType    MeasureType;
    typedef Superclass::DerivativeType DerivativeType;
    typedef vnl_matrix <double> MatrixType;
    typedef Superclass::ParametersType ParametersType;

    /**
     * The measure type shall be used for computing the cost function value to observe convergence
     * Parameters are set as {flip_angle,theta_1,theta_2,theta_3} for unconstrained estimation
     * Or {flip_angle,theta_2} for constrained
     */
    virtual MeasureType GetValue(const ParametersType & parameters) const ITK_OVERRIDE;
    virtual void GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const ITK_OVERRIDE;

    itkSetMacro(EchoSpacing, double)
    itkSetMacro(ExcitationFlipAngle, double)

    void SetT2RelaxometrySignals(ParametersType &relaxoSignals) {m_T2RelaxometrySignals = relaxoSignals;}

    itkSetMacro(T1Value, double)
    void SetGammaMeans(std::vector <double> &val) {m_GammaMeans = val;}
    void SetGammaVariances(std::vector <double> &val) {m_GammaVariances = val;}

    itkSetMacro(GammaIntegralTolerance, double)
    itkSetMacro(ConstrainedParameters, bool)

    unsigned int GetNumberOfParameters() const ITK_OVERRIDE
    {
        if (m_ConstrainedParameters)
            return 2;
        else
            return 4;
    }

    itkGetMacro(SigmaSquare, double)
    ParametersType &GetOptimalT2Weights() {return m_OptimalT2Weights;}

protected:
    B1GammaMixtureT2RelaxometryCostFunction()
    {
            m_NNLSBordersOptimizer = anima::NNLSOptimizer::New();

            m_T1Value = 1;
            m_EchoSpacing = 1;
    }

    virtual ~B1GammaMixtureT2RelaxometryCostFunction() {}

    void PrepareDataForLLS() const;
    void PrepareDataForDerivative() const;

    //! Computes maximum likelihood estimates of weights
    void SolveLinearLeastSquares() const;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(B1GammaMixtureT2RelaxometryCostFunction);

    double m_EchoSpacing;

    ParametersType m_T2RelaxometrySignals;
    mutable ParametersType m_TestedParameters;
    bool m_ConstrainedParameters;

    double m_ExcitationFlipAngle;

    mutable ParametersType m_OptimalT2Weights;
    double m_T1Value;

    mutable std::vector <double> m_GammaMeans, m_GammaVariances;
    double m_GammaIntegralTolerance;

    // Internal working variables, not thread safe but so much faster !
    mutable anima::EPGSignalSimulator m_T2SignalSimulator;

    mutable ParametersType m_FSignals;
    mutable ParametersType m_Residuals;
    mutable vnl_matrix <double> m_PredictedSignalAttenuations, m_CholeskyMatrix;
    mutable double m_SigmaSquare;

    mutable std::vector <MatrixType> m_SignalAttenuationsJacobian;
    mutable MatrixType m_GramMatrix, m_InverseGramMatrix, m_FMatrixInverseG;
    mutable std::vector <bool> m_CompartmentSwitches;

    mutable anima::NNLSOptimizer::Pointer m_NNLSBordersOptimizer;
    mutable anima::CholeskyDecomposition m_CholeskySolver;
};

} // end namespace anima
