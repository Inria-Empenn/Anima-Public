#pragma once

#include <AnimaMCMExport.h>
#include <animaBaseMCMCost.h>

namespace anima {

class ANIMAMCM_EXPORT NonCentralChiMCMCost : public anima::BaseMCMCost {
public:
  /** Standard class typedefs. */
  using Self = NonCentralChiMCMCost;
  using Superclass = anima::BaseMCMCost;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(NonCentralChiMCMCost, anima::BaseMCMCost);

  void SetNumberOfCoils(unsigned int val) { m_NumberOfCoils = val; }

  //! Not really doing anything, only setting parameters for now
  MeasureType GetValues(const ParametersType &parameters) ITK_OVERRIDE;

  //! For the current set of parameters, compute the cost function value,
  //! requires GetValues to be called first
  double GetCurrentCostValue() ITK_OVERRIDE;

  //! Purposedly not implemented
  void GetDerivativeMatrix(const ParametersType &parameters,
                           DerivativeMatrixType &derivative) ITK_OVERRIDE {}

  //! Derivative, not implemented yet
  void GetCurrentDerivative(DerivativeMatrixType &derivativeMatrix,
                            DerivativeType &derivative) ITK_OVERRIDE {
    itkExceptionMacro("Derivative not implemented for non-central chi noise");
  }

protected:
  NonCentralChiMCMCost() { m_NumberOfCoils = 1; }

  virtual ~NonCentralChiMCMCost() {}
  void ComputeSigmaSquareValue();

private:
  NonCentralChiMCMCost(const Self &); // purposely not implemented
  void operator=(const Self &);       // purposely not implemented

  unsigned int m_NumberOfCoils;
};

} // end namespace anima
