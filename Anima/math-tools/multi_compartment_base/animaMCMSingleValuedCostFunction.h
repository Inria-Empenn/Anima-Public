#pragma once

#include <AnimaMCMBaseExport.h>
#include <animaBaseMCMCost.h>
#include <itkSingleValuedCostFunction.h>

namespace anima {

class ANIMAMCMBASE_EXPORT MCMSingleValuedCostFunction
    : public itk::SingleValuedCostFunction {
public:
  /** Standard class typedefs. */
  using Self = MCMSingleValuedCostFunction;
  using Superclass = itk::SingleValuedCostFunction;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MCMSingleValuedCostFunction, itk::SingleValuedCostFunction);

  using MeasureType = Superclass::MeasureType;
  using DerivativeType = Superclass::DerivativeType;
  using ParametersType = Superclass::ParametersType;

  using InternalCostType = anima::BaseMCMCost;
  using InternalCostPointer = InternalCostType::Pointer;

  itkSetMacro(InternalCost, InternalCostPointer);

  virtual MeasureType
  GetValue(const ParametersType &parameters) const ITK_OVERRIDE;
  virtual void GetDerivative(const ParametersType &parameters,
                             DerivativeType &derivative) const ITK_OVERRIDE;

  double GetSigmaSquare();

  InternalCostType::MCMPointer &GetMCMStructure() {
    return m_InternalCost->GetMCMStructure();
  }

  unsigned int GetNumberOfParameters() const ITK_OVERRIDE;

protected:
  MCMSingleValuedCostFunction() { m_InternalCost = ITK_NULLPTR; }

  virtual ~MCMSingleValuedCostFunction() ITK_OVERRIDE {}

private:
  MCMSingleValuedCostFunction(const Self &); // purposely not implemented
  void operator=(const Self &);              // purposely not implemented

  mutable InternalCostPointer m_InternalCost;
};

} // end namespace anima
