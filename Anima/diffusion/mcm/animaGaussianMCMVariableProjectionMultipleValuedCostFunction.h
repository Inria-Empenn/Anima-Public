#pragma once

#include <AnimaMCMExport.h>
#include <animaGaussianMCMVariableProjectionCost.h>
#include <itkMultipleValuedCostFunction.h>

namespace anima {

class ANIMAMCM_EXPORT GaussianMCMVariableProjectionMultipleValuedCostFunction
    : public itk::MultipleValuedCostFunction {
public:
  /** Standard class typedefs. */
  using Self = GaussianMCMVariableProjectionMultipleValuedCostFunction;
  using Superclass = itk::MultipleValuedCostFunction;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(GaussianMCMVariableProjectionMultipleValuedCostFunction,
               itk::MultipleValuedCostFunction);

  using MeasureType = Superclass::MeasureType;
  using DerivativeType = Superclass::DerivativeType;
  using ParametersType = Superclass::ParametersType;

  using InternalCostType = anima::GaussianMCMVariableProjectionCost;
  using InternalCostPointer = InternalCostType::Pointer;

  itkSetMacro(InternalCost, InternalCostPointer);
  itkGetConstReferenceMacro(InternalCost, InternalCostPointer);

  virtual MeasureType
  GetValue(const ParametersType &parameters) const ITK_OVERRIDE;
  virtual void GetDerivative(const ParametersType &parameters,
                             DerivativeType &derivative) const ITK_OVERRIDE;

  double GetSigmaSquare();
  std::vector<double> &GetOptimalWeights();

  unsigned int GetNumberOfParameters() const ITK_OVERRIDE;
  unsigned int GetNumberOfValues() const ITK_OVERRIDE;

protected:
  GaussianMCMVariableProjectionMultipleValuedCostFunction() {
    m_InternalCost = ITK_NULLPTR;
  }

  virtual ~GaussianMCMVariableProjectionMultipleValuedCostFunction() {}

private:
  GaussianMCMVariableProjectionMultipleValuedCostFunction(
      const Self &);            // purposely not implemented
  void operator=(const Self &); // purposely not implemented

  mutable InternalCostPointer m_InternalCost;
};

} // end namespace anima
