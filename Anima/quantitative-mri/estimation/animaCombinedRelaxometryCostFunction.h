#pragma once

#include "AnimaRelaxometryExport.h"
#include <itkSingleValuedCostFunction.h>

namespace anima {

class ANIMARELAXOMETRY_EXPORT CombinedRelaxometryCostFunction
    : public itk::SingleValuedCostFunction {
public:
  /** Standard class typedefs. */
  using Self = CombinedRelaxometryCostFunction;
  using Superclass = itk::SingleValuedCostFunction;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(CombinedRelaxometryCostFunction, Superclass);

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
  itkSetMacro(T2ExcitationFlipAngle, double);

  void SetT1RelaxometrySignals(std::vector<double> &relaxoSignals) {
    m_T1RelaxometrySignals = relaxoSignals;
  }
  void SetT2RelaxometrySignals(std::vector<double> &relaxoSignals) {
    m_T2RelaxometrySignals = relaxoSignals;
  }

  void SetT1FlipAngles(std::vector<double> &flipAngles) {
    m_T1FlipAngles = flipAngles;
  }
  void SetT2FlipAngles(std::vector<double> &flipAngles) {
    m_T2FlipAngles = flipAngles;
  }

  itkSetMacro(TRValue, double);
  itkSetMacro(KValue, double);
  itkSetMacro(T1Value, double);
  itkSetMacro(T2Value, double);
  itkSetMacro(B1Value, double);
  itkSetMacro(B1T2AdditiveValue, double);
  itkSetMacro(M0Value, double);

  enum OptimizedValueType { M0 = 0, B1, B1_Additive, T2, T1 };

  itkSetMacro(OptimizedValue, OptimizedValueType);

  unsigned int GetNumberOfParameters() const ITK_OVERRIDE { return 1; }

protected:
  CombinedRelaxometryCostFunction() {
    m_KValue = 1;
    m_T1Value = 1;
    m_T2Value = 1;
    m_B1Value = 1;
    m_B1T2AdditiveValue = 0;
    m_M0Value = 1;

    m_OptimizedValue = T1;

    m_TRValue = 15;
    m_T2EchoSpacing = 1;
  }

  virtual ~CombinedRelaxometryCostFunction() {}

private:
  CombinedRelaxometryCostFunction(const Self &); // purposely not implemented
  void operator=(const Self &);                  // purposely not implemented

  double m_T2EchoSpacing;
  std::vector<double> m_T1RelaxometrySignals;
  std::vector<double> m_T2RelaxometrySignals;

  double m_T2ExcitationFlipAngle;
  std::vector<double> m_T1FlipAngles;
  std::vector<double> m_T2FlipAngles;

  double m_TRValue;
  double m_KValue;
  mutable double m_T1Value, m_T2Value, m_M0Value, m_B1Value,
      m_B1T2AdditiveValue;

  OptimizedValueType m_OptimizedValue;
};

} // end namespace anima
