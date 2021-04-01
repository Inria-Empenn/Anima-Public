#pragma once
#include <itkSingleValuedCostFunction.h>
#include "AnimaRelaxometryExport.h"

#include <animaEPGSignalSimulator.h>
#include <animaCholeskyDecomposition.h>
#include <animaNNLSOptimizer.h>

#include <map>

namespace anima
{

/**
 * \class B1GMMRelaxometryCostFunction
 * \brief Cost function for estimating B1 from T2 relaxometry acquisition, following a gaussian mixture EPG decay model.
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

    virtual MeasureType GetValue(const ParametersType &parameters) const ITK_OVERRIDE;
    virtual void GetDerivative(const ParametersType &parameters, DerivativeType &derivative) const ITK_OVERRIDE;

    itkSetMacro(EchoSpacing, double)
    itkSetMacro(ExcitationFlipAngle, double)

    void SetT2RelaxometrySignals(ParametersType &relaxoSignals) {m_T2RelaxometrySignals = relaxoSignals;}

    itkSetMacro(T1Value, double)
    void SetGaussianMeans(std::vector <double> &val) {m_GaussianMeans = val;m_TruncatedGaussianIntegrals.clear();}
    void SetGaussianVariances(std::vector <double> &val) {m_GaussianVariances = val;m_TruncatedGaussianIntegrals.clear();}

    itkSetMacro(UniformPulses, bool)
    itkSetMacro(PixelWidth, double)
    void SetPulseProfile(std::vector < std::pair <double, double> > &profile) {m_PulseProfile = profile;}
    void SetExcitationProfile(std::vector < std::pair <double, double> > &profile) {m_ExcitationProfile = profile;}

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

        m_UniformPulses = true;
        m_PixelWidth = 3.0;
    }

    virtual ~B1GMMRelaxometryCostFunction() {}

    void PrepareDataForLLS() const;

    //! Computes maximum likelihood estimates of weights
    void SolveLinearLeastSquares() const;

private:
    B1GMMRelaxometryCostFunction(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    double m_EchoSpacing;
    ParametersType m_T2RelaxometrySignals;
    mutable ParametersType m_TestedParameters;

    double m_ExcitationFlipAngle;

    mutable ParametersType m_OptimalT2Weights;
    double m_T1Value;
    std::vector <double> m_GaussianMeans, m_GaussianVariances;
    mutable std::vector <double> m_TruncatedGaussianIntegrals;

    // Internal working variables, not thread safe but so much faster !
    mutable anima::EPGSignalSimulator m_T2SignalSimulator;

    bool m_UniformPulses;
    std::vector < std::pair <double, double> > m_PulseProfile;
    std::vector < std::pair <double, double> > m_ExcitationProfile;
    double m_PixelWidth;

    mutable ParametersType m_FSignals;
    mutable ParametersType m_Residuals;
    mutable vnl_matrix <double> m_PredictedSignalAttenuations, m_CholeskyMatrix;
    mutable double m_SigmaSquare;

    mutable anima::NNLSOptimizer::Pointer m_NNLSBordersOptimizer;
    mutable anima::CholeskyDecomposition m_CholeskySolver;
};

} // end namespace anima
