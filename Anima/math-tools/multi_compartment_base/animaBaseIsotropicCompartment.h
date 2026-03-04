#pragma once

#include <AnimaMCMBaseExport.h>
#include <animaBaseCompartment.h>

namespace anima {

class ANIMAMCMBASE_EXPORT BaseIsotropicCompartment : public BaseCompartment {
public:
  // Useful typedefs
  using Self = BaseIsotropicCompartment;
  using Superclass = BaseCompartment;
  using BasePointer = Superclass::Pointer;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;
  using ModelOutputVectorType = Superclass::ModelOutputVectorType;
  using Vector3DType = Superclass::Vector3DType;
  using Matrix3DType = Superclass::Matrix3DType;

  /** Run-time type information (and related methods) */
  itkTypeMacro(BaseIsotropicCompartment, BaseCompartment);

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

  // Set constraints
  void SetEstimateAxialDiffusivity(bool arg);
  bool GetEstimateAxialDiffusivity() { return m_EstimateAxialDiffusivity; }

  void
  SetCompartmentVector(ModelOutputVectorType &compartmentVector) ITK_OVERRIDE;

  unsigned int GetCompartmentSize() ITK_OVERRIDE;
  unsigned int GetNumberOfParameters() ITK_OVERRIDE;
  ModelOutputVectorType &GetCompartmentVector() ITK_OVERRIDE;

  //! Reimplements re-orientation, useless in isotropic water compartment
  void Reorient(vnl_matrix<double> &orientationMatrix,
                bool affineTransform) ITK_OVERRIDE {}

  const Matrix3DType &GetDiffusionTensor() ITK_OVERRIDE;
  double GetApparentFractionalAnisotropy() ITK_OVERRIDE;
  double GetApparentMeanDiffusivity() ITK_OVERRIDE;
  double GetApparentParallelDiffusivity() ITK_OVERRIDE;
  double GetApparentPerpendicularDiffusivity() ITK_OVERRIDE;

protected:
  BaseIsotropicCompartment() : Superclass() {
    m_EstimateAxialDiffusivity = true;
    m_ChangedConstraints = true;

    m_NumberOfParameters = this->GetCompartmentSize();
  }

  virtual ~BaseIsotropicCompartment() {}

private:
  bool m_EstimateAxialDiffusivity;
  bool m_ChangedConstraints;
  unsigned int m_NumberOfParameters;

  Matrix3DType m_DiffusionTensor;
};

} // end namespace anima
