#pragma once

#include <itkSingleValuedCostFunction.h>
#include "AnimaRelaxometryExport.h"

#include <animaEPGSignalSimulator.h>

namespace anima
{
    
/** \class B1T2RelaxometryDistributionCostFunction
 * \brief Cost function for estimating B1 from T2 relaxometry acquisition, following a multi-T2 EPG decay model.
 *
 */
class ANIMARELAXOMETRY_EXPORT B1T2RelaxometryDistributionCostFunction :
public itk::SingleValuedCostFunction
{
public:
    /** Standard class typedefs. */
    typedef B1T2RelaxometryDistributionCostFunction Self;
    typedef itk::SingleValuedCostFunction   Superclass;
    typedef itk::SmartPointer<Self>         Pointer;
    typedef itk::SmartPointer<const Self>   ConstPointer;

    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(B1T2RelaxometryDistributionCostFunction, Superclass)

    typedef Superclass::MeasureType MeasureType;
    typedef Superclass::DerivativeType DerivativeType;
    typedef Superclass::ParametersType ParametersType;
    typedef std::vector < std::complex <double> > ComplexVectorType;
    typedef std::vector <ComplexVectorType> MatrixType;

    virtual MeasureType GetValue(const ParametersType & parameters) const ITK_OVERRIDE;
    virtual void GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const ITK_OVERRIDE {} // No derivative

    itkSetMacro(EchoSpacing, double)
    itkSetMacro(ExcitationFlipAngle, double)

    void SetT2RelaxometrySignals(ParametersType &relaxoSignals) {m_T2RelaxometrySignals = relaxoSignals;}
    void SetT2FlipAngles(std::vector <double> & flipAngles) {m_T2FlipAngles = flipAngles;}

    itkSetMacro(T1Value, double)
    itkSetMacro(M0Value, double)

    void SetT2DistributionSamples(std::vector < std::vector <double> > &values) {m_T2DistributionSamples = values;}
    void SetLowerT2Bound(double value) {m_LowerT2Bound = value;}
    void SetUpperT2Bound(double value) {m_UpperT2Bound = value;}
    void SetT2IntegrationStep(double value) {m_T2IntegrationStep = value;}

    void SetT2Weights(ParametersType &weights) {m_T2Weights = weights;}

    itkSetMacro(B1OnExcitationAngle, bool)

    void SetT2WorkingValues(std::vector <double> &values) {m_T2WorkingValues = values;}
    void SetDistributionSamplesT2Correspondences(std::vector < std::vector <unsigned int> > &values) {m_DistributionSamplesT2Correspondences = values;}

    unsigned int GetNumberOfParameters() const ITK_OVERRIDE
    {
        return 1;
    }

protected:
    B1T2RelaxometryDistributionCostFunction()
    {
        m_T1Value = 1;
        m_M0Value = 1;

        m_EchoSpacing = 1;

        m_T2IntegrationStep = 1;
        m_LowerT2Bound = 1.0e-4;
        m_UpperT2Bound = 3000;

        m_B1OnExcitationAngle = false;
    }

    virtual ~B1T2RelaxometryDistributionCostFunction() {}

private:
    B1T2RelaxometryDistributionCostFunction(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    double m_EchoSpacing;
    ParametersType m_T2RelaxometrySignals;

    double m_ExcitationFlipAngle;
    std::vector <double> m_T2FlipAngles;

    double m_T2IntegrationStep;
    double m_LowerT2Bound, m_UpperT2Bound;
    std::vector < std::vector <double> > m_T2DistributionSamples;
    std::vector <double> m_T2WorkingValues;
    std::vector < std::vector <unsigned int> > m_DistributionSamplesT2Correspondences;

    ParametersType m_T2Weights;

    bool m_B1OnExcitationAngle;
    double m_T1Value, m_M0Value;

    // Internal working variables, not thread safe but so much faster !
    mutable anima::EPGSignalSimulator m_T2SignalSimulator;
    mutable std::vector <anima::EPGSignalSimulator::RealVectorType> m_SimulatedEPGValues;
    mutable anima::EPGSignalSimulator::RealVectorType m_SimulatedSignalValues;
};
    
} // end namespace anima
