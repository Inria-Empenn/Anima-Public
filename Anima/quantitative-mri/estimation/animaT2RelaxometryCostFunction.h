#pragma once

#include "AnimaRelaxometryExport.h"
#include <itkSingleValuedCostFunction.h>

namespace anima {

class ANIMARELAXOMETRY_EXPORT T2RelaxometryCostFunction
    : public itk::SingleValuedCostFunction {
public:
  /** Standard class typedefs. */
  using Self = T2RelaxometryCostFunction;
  using Superclass = itk::SingleValuedCostFunction;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(T2RelaxometryCostFunction, Superclass);

  using MeasureType = Superclass::MeasureType;
  using DerivativeType = Superclass::DerivativeType;
  using ParametersType = Superclass::ParametersType;
  using ComplexVectorType = std::vector<std::complex<double>>;
  using MatrixType = std::vector<ComplexVectorType>;

  virtual MeasureType
  GetValue(const ParametersType &parameters) const ITK_OVERRIDE;
  virtual void GetDerivative(const ParametersType &parameters,
                             DerivativeType &derivative) const ITK_OVERRIDE {
  } // No derivative

  itkSetMacro(T2EchoSpacing, double);

  void SetT2RelaxometrySignals(std::vector<double> &relaxoSignals) {
    m_T2RelaxometrySignals = relaxoSignals;
  }

  itkSetMacro(T1Value, double);
  itkSetMacro(T2Value, double);
  itkSetMacro(TRValue, double);
  itkGetMacro(M0Value, double);

  unsigned int GetNumberOfParameters() const ITK_OVERRIDE {
    // T2
    return 1;
  }

protected:
  T2RelaxometryCostFunction() {
    m_T1Value = 1;
    m_T2Value = 1;
    m_M0Value = 1;
    m_TRValue = 5000;

    m_T2EchoSpacing = 1;
  }

  virtual ~T2RelaxometryCostFunction() {}

private:
  T2RelaxometryCostFunction(const Self &); // purposely not implemented
  void operator=(const Self &);            // purposely not implemented

  double m_T2EchoSpacing;
  std::vector<double> m_T2RelaxometrySignals;

  double m_TRValue;
  mutable double m_T1Value, m_T2Value, m_M0Value;
};

} // end namespace anima
