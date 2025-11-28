#pragma once

#include <itkArray2D.h>
#include <itkLightObject.h>
#include <itkOptimizerParameters.h>

#include <AnimaMCMBaseExport.h>
#include <animaMultiCompartmentModel.h>

namespace anima {

/**
 * @brief Base cost function class to handle maximum likelihood estimation
 */
class ANIMAMCMBASE_EXPORT BaseMCMCost : public itk::LightObject {
public:
  /** Standard class typedefs. */
  using Self = BaseMCMCost;
  using Superclass = itk::LightObject;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Run-time type information (and related methods). */
  itkTypeMacro(BaseMCMCost, itk::LightObject);

  using DerivativeMatrixType = itk::Array2D<double>;
  using MeasureType = itk::Array<double>;
  using DerivativeType = itk::Array<double>;
  using ParametersType = itk::OptimizerParameters<double>;

  using MCMType = anima::MultiCompartmentModel;
  using MCMPointer = MCMType::Pointer;
  using Vector3DType = MCMType::Vector3DType;
  using ListType = MCMType::ListType;

  void SetObservedSignals(ListType &value) { m_ObservedSignals = value; }
  void SetGradients(std::vector<Vector3DType> &value) { m_Gradients = value; }
  void SetGradientStrengths(ListType &value) { m_GradientStrengths = value; }

  void SetMCMStructure(MCMType *model) { m_MCMStructure = model; }
  MCMPointer &GetMCMStructure() { return m_MCMStructure; }

  //! Get residual values for a given set of parameters, returns a vector of
  //! residuals
  virtual MeasureType GetValues(const ParametersType &parameters) = 0;

  //! For the current set of parameters, compute the cost function value,
  //! requires GetValues to be called first
  virtual double GetCurrentCostValue() = 0;

  //! Get residual derivatives for a given set of parameters, returns a matrix
  //! of residuals derivatives
  virtual void GetDerivativeMatrix(const ParametersType &parameters,
                                   DerivativeMatrixType &derivative) = 0;

  //! Get cost function derivatives from a derivative matrix obtained from
  //! GetDerivativeMatrix
  virtual void GetCurrentDerivative(DerivativeMatrixType &derivativeMatrix,
                                    DerivativeType &derivative) = 0;

  //! Returns number of optimized parameters
  unsigned int GetNumberOfParameters() const {
    if (!m_MCMStructure)
      itkExceptionMacro("A multi-compartment model has to be set, no "
                        "computation can be done otherwise.");

    return m_MCMStructure->GetNumberOfParameters();
  }

  unsigned int GetNumberOfObservations() const {
    return m_ObservedSignals.size();
  }

  virtual double GetSigmaSquare() { return m_SigmaSquare; }

  void SetSmallDelta(double val) { m_SmallDelta = val; }
  void SetBigDelta(double val) { m_BigDelta = val; }

protected:
  BaseMCMCost();
  virtual ~BaseMCMCost() {}

  double m_SigmaSquare;
  std::vector<double> m_PredictedSignals;

  ListType m_TestedParameters;

  ListType m_ObservedSignals;
  std::vector<Vector3DType> m_Gradients;

  double m_SmallDelta;
  double m_BigDelta;
  ListType m_GradientStrengths;

  MCMPointer m_MCMStructure;

private:
  BaseMCMCost(const Self &);    // purposely not implemented
  void operator=(const Self &); // purposely not implemented
};

} // end namespace anima
