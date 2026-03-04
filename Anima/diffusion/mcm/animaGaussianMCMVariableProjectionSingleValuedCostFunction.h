#pragma once

#include <AnimaMCMExport.h>
#include <animaGaussianMCMVariableProjectionCost.h>
#include <itkSingleValuedCostFunction.h>

namespace anima {

class ANIMAMCM_EXPORT GaussianMCMVariableProjectionSingleValuedCostFunction
    : public itk::SingleValuedCostFunction {
public:
  /** Standard class typedefs. */
  using Self = GaussianMCMVariableProjectionSingleValuedCostFunction;
  using Superclass = itk::SingleValuedCostFunction;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(GaussianMCMVariableProjectionSingleValuedCostFunction,
               itk::SingleValuedCostFunction);

  using MeasureType = Superclass::MeasureType;
  using DerivativeType = Superclass::DerivativeType;
  using ParametersType = Superclass::ParametersType;

  using InternalCostType = anima::GaussianMCMVariableProjectionCost;
  using InternalCostPointer = InternalCostType::Pointer;

  itkSetMacro(InternalCost, InternalCostPointer);

  virtual MeasureType
  GetValue(const ParametersType &parameters) const ITK_OVERRIDE;
  virtual void GetDerivative(const ParametersType &parameters,
                             DerivativeType &derivative) const ITK_OVERRIDE;

  double GetSigmaSquare();
  std::vector<double> &GetOptimalWeights();

  unsigned int GetNumberOfParameters() const ITK_OVERRIDE;

protected:
  GaussianMCMVariableProjectionSingleValuedCostFunction() {
    m_InternalCost = ITK_NULLPTR;
  }

  virtual ~GaussianMCMVariableProjectionSingleValuedCostFunction()
      ITK_OVERRIDE {}

private:
  GaussianMCMVariableProjectionSingleValuedCostFunction(
      const Self &);            // purposely not implemented
  void operator=(const Self &); // purposely not implemented

  mutable InternalCostPointer m_InternalCost;
};

} // end namespace anima
