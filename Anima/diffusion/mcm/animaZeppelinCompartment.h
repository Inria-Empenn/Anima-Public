#pragma once

#include <AnimaMCMExport.h>
#include <animaBaseCompartment.h>
#include <animaMatrixOperations.h>

namespace anima {

class ANIMAMCM_EXPORT ZeppelinCompartment : public BaseCompartment {
public:
  // Useful typedefs
  using Self = ZeppelinCompartment;
  using Superclass = BaseCompartment;
  using BasePointer = Superclass::Pointer;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;
  using ModelOutputVectorType = Superclass::ModelOutputVectorType;
  using Vector3DType = Superclass::Vector3DType;
  using Matrix3DType = Superclass::Matrix3DType;

  // New macro
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(ZeppelinCompartment, BaseCompartment);

  DiffusionModelCompartmentType GetCompartmentType() ITK_OVERRIDE {
    return Zeppelin;
  }

  virtual double GetFourierTransformedDiffusionProfile(
      double smallDelta, double bigDelta, double gradientStrength,
      const Vector3DType &gradient) ITK_OVERRIDE;
  virtual ListType &
  GetSignalAttenuationJacobian(double smallDelta, double bigDelta,
                               double gradientStrength,
                               const Vector3DType &gradient) ITK_OVERRIDE;
  virtual double
  GetLogDiffusionProfile(const Vector3DType &sample) ITK_OVERRIDE;

  virtual void SetParametersFromVector(const ListType &params) ITK_OVERRIDE;
  virtual ListType &GetParametersAsVector() ITK_OVERRIDE;

  virtual ListType &GetParameterLowerBounds() ITK_OVERRIDE;
  virtual ListType &GetParameterUpperBounds() ITK_OVERRIDE;

  // Set constraints
  void SetEstimateDiffusivities(bool arg);
  void
  SetCompartmentVector(ModelOutputVectorType &compartmentVector) ITK_OVERRIDE;

  unsigned int GetCompartmentSize() ITK_OVERRIDE;
  unsigned int GetNumberOfParameters() ITK_OVERRIDE;
  ModelOutputVectorType &GetCompartmentVector() ITK_OVERRIDE;

  void SetRadialDiffusivity1(double num) ITK_OVERRIDE;

  const Matrix3DType &GetDiffusionTensor() ITK_OVERRIDE;
  double GetApparentFractionalAnisotropy() ITK_OVERRIDE;
  double GetApparentMeanDiffusivity() ITK_OVERRIDE;
  double GetApparentParallelDiffusivity() ITK_OVERRIDE;
  double GetApparentPerpendicularDiffusivity() ITK_OVERRIDE;

protected:
  ZeppelinCompartment() : Superclass() {
    m_EstimateDiffusivities = true;
    m_ChangedConstraints = true;
    m_GradientEigenvector1 = 0;
  }

  virtual ~ZeppelinCompartment() {}

private:
  bool m_EstimateDiffusivities;
  bool m_ChangedConstraints;
  unsigned int m_NumberOfParameters;
  double m_GradientEigenvector1;
};

} // end namespace anima
