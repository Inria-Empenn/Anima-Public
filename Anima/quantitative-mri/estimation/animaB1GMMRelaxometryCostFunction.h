#pragma once

#include <itkSingleValuedCostFunction.h>
#include "AnimaRelaxometryExport.h"
#include <vnl_matrix.h>

#include <animaEPGSignalSimulator.h>
#include <animaCholeskyDecomposition.h>
#include <animaNNLSOptimizer.h>

namespace anima
{
    
/** \class B1GMMRelaxometryCostFunction
 * \brief Cost function for estimating B1 from T2 relaxometry acquisition, following a multi-T2 EPG decay model.
 * The cost function includes (via variable projection) estimation of compartment weights
 *
 */
class ANIMARELAXOMETRY_EXPORT B1GMMRelaxometryCostFunction :
public itk::SingleValuedCostFunction
{
public:
    /** Standard class typedefs. */
    typedef B1GMMRelaxometryCostFunction Self;
    typedef itk::SingleValuedCostFunction Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(B1GMMRelaxometryCostFunction, Superclass)

    typedef Superclass::MeasureType MeasureType;
    typedef Superclass::DerivativeType DerivativeType;
    typedef Superclass::ParametersType ParametersType;
    typedef std::vector < std::complex <double> > ComplexVectorType;
    typedef std::vector <ComplexVectorType> MatrixType;

    virtual MeasureType GetValue(const ParametersType & parameters) const ITK_OVERRIDE;
    virtual void GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const ITK_OVERRIDE;

    itkSetMacro(EchoSpacing, double)
    itkSetMacro(UseDerivative, bool)
    itkSetMacro(ExcitationFlipAngle, double)

    void SetT2RelaxometrySignals(ParametersType &relaxoSignals) {m_T2RelaxometrySignals = relaxoSignals;}
    void SetT2FlipAngles(std::vector <double> & flipAngles) {m_T2FlipAngles = flipAngles;}

    itkSetMacro(T1Value, double)

    void SetT2DistributionSamples(std::vector < std::vector <double> > &values) {m_T2DistributionSamples = values;}
    void SetLowerT2Bound(double value) {m_LowerT2Bound = value;}
    void SetUpperT2Bound(double value) {m_UpperT2Bound = value;}
    void SetT2IntegrationStep(double value) {m_T2IntegrationStep = value;}

    void SetT2WorkingValues(std::vector <double> &values) {m_T2WorkingValues = values;}
    void SetDistributionSamplesT2Correspondences(std::vector < std::vector <unsigned int> > &values) {m_DistributionSamplesT2Correspondences = values;}

    unsigned int GetNumberOfParameters() const ITK_OVERRIDE
    {
        return 1;
    }

    itkGetMacro(SigmaSquare, double)
    ParametersType &GetOptimalT2Weights() {return m_OptimalT2Weights;}

protected:
    B1GMMRelaxometryCostFunction()
    {
        m_NNLSBordersOptimizer = anima::NNLSOptimizer::New();

        m_T1Value = 1;
        m_EchoSpacing = 1;
        m_UseDerivative = true;

        m_T2IntegrationStep = 1;
        m_LowerT2Bound = 1.0e-4;
        m_UpperT2Bound = 3000;
    }

    virtual ~B1GMMRelaxometryCostFunction() {}

    void PrepareDataForLLS() const;
    void PrepareDataForDerivative() const;

    //! Computes maximum likelihood estimates of weights
    void SolveLinearLeastSquares() const;

    void GetDerivativeMatrix(const ParametersType &parameters, DerivativeType &derivative) const;

private:
    B1GMMRelaxometryCostFunction(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    double m_EchoSpacing;
    bool m_UseDerivative;
    ParametersType m_T2RelaxometrySignals;
    mutable ParametersType m_TestedParameters;

    double m_ExcitationFlipAngle;
    std::vector <double> m_T2FlipAngles;

    double m_T2IntegrationStep;
    double m_LowerT2Bound, m_UpperT2Bound;
    std::vector < std::vector <double> > m_T2DistributionSamples;
    std::vector <double> m_T2WorkingValues;
    std::vector < std::vector <unsigned int> > m_DistributionSamplesT2Correspondences;

    mutable ParametersType m_OptimalT2Weights;
    double m_T1Value;

    // Internal working variables, not thread safe but so much faster !
    mutable anima::EPGSignalSimulator m_T2SignalSimulator;
    mutable std::vector <anima::EPGSignalSimulator::RealVectorType> m_SimulatedEPGValues;
    mutable std::vector <anima::EPGSignalSimulator::RealVectorType> m_SimulatedEPGDerivatives;
    mutable anima::EPGSignalSimulator::RealVectorType m_SimulatedSignalValues;

    mutable ParametersType m_FSignals;
    mutable vnl_matrix <double> m_FMatrixInverseG;
    mutable std::vector <bool> m_CompartmentSwitches;
    mutable ParametersType m_Residuals;
    mutable vnl_matrix <double> m_SignalAttenuationsJacobian;
    mutable vnl_matrix <double> m_PredictedSignalAttenuations, m_CholeskyMatrix;
    mutable double m_SigmaSquare;

    mutable anima::NNLSOptimizer::Pointer m_NNLSBordersOptimizer;
    mutable anima::CholeskyDecomposition m_CholeskySolver;
    mutable vnl_matrix <double> m_GramMatrix, m_InverseGramMatrix;
};
    
} // end namespace anima
