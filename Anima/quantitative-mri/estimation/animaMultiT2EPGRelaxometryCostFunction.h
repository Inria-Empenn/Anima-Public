#pragma once

#include <vnl/vnl_matrix.h>
#include <itkSingleValuedCostFunction.h>
#include <animaNNLSOptimizer.h>
#include "AnimaRelaxometryExport.h"

namespace anima
{
    
/** \class MultiT2EPGRelaxometryCostFunction
 * \brief Cost function for estimating B1 from T2 relaxometry acquisition, following a multi-T2 EPG decay model.
 *
 */
class ANIMARELAXOMETRY_EXPORT MultiT2EPGRelaxometryCostFunction :
public itk::SingleValuedCostFunction
{
public:
    /** Standard class typedefs. */
    typedef MultiT2EPGRelaxometryCostFunction Self;
    typedef itk::SingleValuedCostFunction   Superclass;
    typedef itk::SmartPointer<Self>         Pointer;
    typedef itk::SmartPointer<const Self>   ConstPointer;

    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(MultiT2EPGRelaxometryCostFunction, Superclass);

    typedef Superclass::MeasureType MeasureType;
    typedef Superclass::DerivativeType DerivativeType;
    typedef Superclass::ParametersType ParametersType;
    typedef std::vector < std::complex <double> > ComplexVectorType;
    typedef std::vector <ComplexVectorType> MatrixType;

    typedef anima::NNLSOptimizer NNLSOptimizerType;
    typedef NNLSOptimizerType::Pointer NNLSOptimizerPointer;

    virtual MeasureType GetValue(const ParametersType & parameters) const ITK_OVERRIDE;
    virtual void GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const ITK_OVERRIDE {} // No derivative

    itkSetMacro(EchoSpacing, double)
    itkSetMacro(ExcitationFlipAngle, double)

    void SetT2RelaxometrySignals(ParametersType &relaxoSignals) {m_T2RelaxometrySignals = relaxoSignals;}

    itkSetMacro(T1Value, double)
    itkGetMacro(OptimizedM0Value, double)

    void SetT2Values(std::vector <double> &values) {m_T2Values = values;}
    ParametersType &GetOptimizedT2Weights() {return m_OptimizedT2Weights;}
    vnl_matrix <double> &GetAMatrix() {return m_AMatrix;}

    unsigned int GetNumberOfParameters() const ITK_OVERRIDE
    {
        return 1;
    }

protected:
    MultiT2EPGRelaxometryCostFunction()
    {
        m_T1Value = 1;
        m_OptimizedM0Value = 1;
        m_EchoSpacing = 1;

        m_NNLSOptimizer = NNLSOptimizerType::New();
    }

    virtual ~MultiT2EPGRelaxometryCostFunction() {}

private:
    MultiT2EPGRelaxometryCostFunction(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    double m_EchoSpacing;
    ParametersType m_T2RelaxometrySignals;

    double m_ExcitationFlipAngle;
    std::vector <double> m_T2Values;

    double m_T1Value;

    mutable NNLSOptimizerPointer m_NNLSOptimizer;
    mutable vnl_matrix <double> m_AMatrix;
    mutable ParametersType m_OptimizedT2Weights;
    mutable double m_OptimizedM0Value;
};
    
} // end namespace anima
