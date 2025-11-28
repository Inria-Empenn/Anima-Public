#pragma once

#include <AnimaMCMBaseExport.h>
#include <animaBaseMCMCost.h>
#include <itkMultipleValuedCostFunction.h>

namespace anima {

class ANIMAMCMBASE_EXPORT MCMMultipleValuedCostFunction
    : public itk::MultipleValuedCostFunction {
public:
  /** Standard class typedefs. */
  using Self = MCMMultipleValuedCostFunction;
  using Superclass = itk::MultipleValuedCostFunction;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MCMMultipleValuedCostFunction, itk::MultipleValuedCostFunction);

  using MeasureType = Superclass::MeasureType;
  using DerivativeType = Superclass::DerivativeType;
  using ParametersType = Superclass::ParametersType;

  using InternalCostType = anima::BaseMCMCost;
  using InternalCostPointer = InternalCostType::Pointer;

  itkSetMacro(InternalCost, InternalCostPointer);
  itkGetConstReferenceMacro(InternalCost, InternalCostPointer);

  virtual MeasureType
  GetValue(const ParametersType &parameters) const ITK_OVERRIDE;
  virtual void GetDerivative(const ParametersType &parameters,
                             DerivativeType &derivative) const ITK_OVERRIDE;

  double GetSigmaSquare();

  InternalCostType::MCMPointer &GetMCMStructure() {
    return m_InternalCost->GetMCMStructure();
  }

  unsigned int GetNumberOfParameters() const ITK_OVERRIDE;
  unsigned int GetNumberOfValues() const ITK_OVERRIDE;

protected:
  MCMMultipleValuedCostFunction() { m_InternalCost = ITK_NULLPTR; }

  virtual ~MCMMultipleValuedCostFunction() ITK_OVERRIDE {}

private:
  MCMMultipleValuedCostFunction(const Self &); // purposely not implemented
  void operator=(const Self &);                // purposely not implemented

  mutable InternalCostPointer m_InternalCost;
};

} // end namespace anima
