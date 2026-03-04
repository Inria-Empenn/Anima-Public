#pragma once

#include <animaMultiCompartmentModel.h>
#include <itkLightObject.h>
#include <itkVariableLengthVector.h>

#include <vnl/vnl_math.h>

#include <AnimaMCMBaseExport.h>

namespace anima {

/**
 * @brief Computes a L2 distance between two MCM of any type
 */
class ANIMAMCMBASE_EXPORT MCML2DistanceComputer : public itk::LightObject {
public:
  using Self = MCML2DistanceComputer;
  using Superclass = itk::LightObject;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Run-time type information (and related methods) */
  itkTypeMacro(MCML2DistanceComputer, itk::LightObject);

  itkNewMacro(Self);

  using MCMType = anima::MultiCompartmentModel;
  using GradientType = MCMType::Vector3DType;
  using MCMPointer = MCMType::Pointer;

  void SetLowPassGaussianSigma(double val) { m_LowPassGaussianSigma = val; }
  void SetForceApproximation(bool val) { m_ForceApproximation = val; }
  void SetSquaredDistance(bool val) { m_SquaredDistance = val; }

  double ComputeDistance(const MCMPointer &firstModel,
                         const MCMPointer &secondModel) const;

  void SetGradientStrengths(const std::vector<double> &val);
  void SetSmallDelta(double val);
  void SetBigDelta(double val);

  void SetGradientDirections(const std::vector<GradientType> &val);

protected:
  MCML2DistanceComputer();
  ~MCML2DistanceComputer() {}

  void UpdateSphereWeights();
  bool CheckTensorCompatibility(const MCMPointer &firstModel,
                                const MCMPointer &secondModel) const;

  double ComputeTensorDistance(const MCMPointer &firstModel,
                               const MCMPointer &secondModel) const;
  double ComputeApproximateDistance(const MCMPointer &firstModel,
                                    const MCMPointer &secondModel) const;

private:
  double m_LowPassGaussianSigma;
  bool m_ForceApproximation;
  bool m_SquaredDistance;

  // Optional parameters for the case when compartments are not tensor
  // compatible
  double m_SmallDelta;
  double m_BigDelta;
  std::vector<double> m_GradientStrengths;
  std::vector<GradientType> m_GradientDirections;

  // Parameters for numerical integration on non tensor compatible models
  std::vector<double> m_SphereWeights;
  std::vector<unsigned int> m_BValWeightsIndexes;
};

} // end namespace anima
