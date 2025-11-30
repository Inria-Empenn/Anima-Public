#pragma once

#include "AnimaRelaxometryExport.h"
#include <itkSingleValuedCostFunction.h>

namespace anima {

class ANIMARELAXOMETRY_EXPORT T1SERelaxometryCostFunction
    : public itk::SingleValuedCostFunction {
public:
  /** Standard class typedefs. */
  using Self = T1SERelaxometryCostFunction;
  using Superclass = itk::SingleValuedCostFunction;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(T1SERelaxometryCostFunction, Superclass);

  using MeasureType = Superclass::MeasureType;
  using DerivativeType = Superclass::DerivativeType;
  using ParametersType = Superclass::ParametersType;
  using ComplexVectorType = std::vector<std::complex<double>>;
  using MatrixType = std::vector<ComplexVectorType>;

  virtual MeasureType
  GetValue(const ParametersType &parameters) const ITK_OVERRIDE;
  virtual void GetDerivative(const ParametersType &parameters,
                             DerivativeType &derivative) const ITK_OVERRIDE;

  void SetRelaxometrySignals(std::vector<double> &relaxoSignals) {
    m_RelaxometrySignals = relaxoSignals;
  }
  void SetTRValues(std::vector<double> &trValues) { m_TRValues = trValues; }

  unsigned int GetNumberOfParameters() const ITK_OVERRIDE { return 2; }

protected:
  T1SERelaxometryCostFunction() {}

  virtual ~T1SERelaxometryCostFunction() {}

private:
  T1SERelaxometryCostFunction(const Self &); // purposely not implemented
  void operator=(const Self &);              // purposely not implemented

  std::vector<double> m_RelaxometrySignals;
  std::vector<double> m_TRValues;
};

} // end namespace anima
