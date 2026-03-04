#pragma once

#include <AnimaMCMSimilarityExport.h>
#include <animaMCMConstants.h>
#include <animaMultiCompartmentModel.h>
#include <itkSingleValuedCostFunction.h>

namespace anima {

class ANIMAMCMSIMILARITY_EXPORT ApproximateMCMSmoothingCostFunction
    : public itk::SingleValuedCostFunction {
public:
  /** Standard class typedefs. */
  using Self = ApproximateMCMSmoothingCostFunction;
  using Superclass = itk::SingleValuedCostFunction;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ApproximateMCMSmoothingCostFunction,
               itk::SingleValuedCostFunction);

  using MeasureType = Superclass::MeasureType;
  using DerivativeType = Superclass::DerivativeType;
  using ParametersType = Superclass::ParametersType;

  using MCModelType = anima::MultiCompartmentModel;
  using MCMPointer = MCModelType::Pointer;
  using GradientType = MCModelType::Vector3DType;
  using BaseCompartmentPointer = anima::BaseCompartment::Pointer;
  using Vector3DType = anima::BaseCompartment::Vector3DType;

  virtual MeasureType
  GetValue(const ParametersType &parameters) const ITK_OVERRIDE;
  virtual void GetDerivative(const ParametersType &parameters,
                             DerivativeType &derivative) const ITK_OVERRIDE;

  void SetReferenceModels(const std::vector<MCMPointer> &refModels,
                          const std::vector<GradientType> &gradients,
                          const double &smallDelta, const double &bigDelta,
                          const std::vector<double> &gradientStrengths);
  void SetMovingModels(const std::vector<MCMPointer> &movingModels,
                       const std::vector<GradientType> &gradients,
                       const double &smallDelta, const double &bigDelta,
                       const std::vector<double> &gradientStrengths);

  void SetGradientStrengths(const std::vector<double> &val);
  void SetGradientDirections(const std::vector<GradientType> &val);

  void SetBValueWeightIndexes(const std::vector<unsigned int> &val);
  void SetSphereWeights(const std::vector<double> &val);

  itkSetMacro(ParameterScale, double);
  void SetSmallDelta(double val);
  void SetBigDelta(double val);

  unsigned int GetNumberOfParameters() const ITK_OVERRIDE { return 1; }

protected:
  ApproximateMCMSmoothingCostFunction() {
    m_UpdatedData = false;
    m_ConstantTerm = 0;
    m_ParameterScale = 1.0e-3;
    m_SmallDelta = anima::DiffusionSmallDelta;
    m_BigDelta = anima::DiffusionBigDelta;
  }

  virtual ~ApproximateMCMSmoothingCostFunction() {}

private:
  ApproximateMCMSmoothingCostFunction(const Self &); // purposely not
                                                     // implemented
  void operator=(const Self &); // purposely not implemented

  std::vector<std::vector<double>> m_ReferenceModelSignalValues;
  std::vector<std::vector<double>> m_MovingModelSignalValues;

  std::vector<unsigned int> m_BValueWeightIndexes;
  std::vector<double> m_GradientStrengths, m_SphereWeights;
  std::vector<GradientType> m_GradientDirections;

  mutable bool m_UpdatedData;
  mutable double m_ConstantTerm;
  double m_ParameterScale;
  double m_SmallDelta, m_BigDelta;
};

} // end namespace anima
