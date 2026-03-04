#pragma once

#include <iostream>
#include <itkMatrixOffsetTransformBase.h>

namespace anima {

template <class TScalarType =
              double> // Data type for scalars (double or double)
class SymmetryPlaneTransform
    : public itk::MatrixOffsetTransformBase<TScalarType, 3> {
public:
  /** Standard class typedefs. */
  using Self = SymmetryPlaneTransform;
  using Superclass = itk::MatrixOffsetTransformBase<TScalarType, 3>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** New macro for creation of through a Smart Pointer. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SymmetryPlaneTransform, MatrixOffsetTransformBase);

  /** Dimension of the space. */
  itkStaticConstMacro(SpaceDimension, unsigned int, 3);
  itkStaticConstMacro(InputSpaceDimension, unsigned int, 3);
  itkStaticConstMacro(OutputSpaceDimension, unsigned int, 3);
  itkStaticConstMacro(ParametersDimension, unsigned int, 3);

  // Check which typedefs are useful
  using ParametersType = typename Superclass::ParametersType;
  using ParametersValueType = typename Superclass::ParametersValueType;
  using JacobianType = typename Superclass::JacobianType;
  using ScalarType = typename Superclass::ScalarType;
  using InputVectorType = typename Superclass::InputVectorType;
  using OutputVectorType = typename Superclass::OutputVectorType;
  using InputCovariantVectorType =
      typename Superclass::InputCovariantVectorType;
  using OutputCovariantVectorType =
      typename Superclass::OutputCovariantVectorType;
  using InputVnlVectorType = typename Superclass::InputVnlVectorType;
  using OutputVnlVectorType = typename Superclass::OutputVnlVectorType;
  using InputPointType = typename Superclass::InputPointType;
  using OutputPointType = typename Superclass::OutputPointType;
  using MatrixType = typename Superclass::MatrixType;
  using HomogeneousMatrixType = itk::Matrix<ScalarType, 4, 4>;
  using InverseMatrixType = typename Superclass::InverseMatrixType;
  using CenterType = typename Superclass::CenterType;
  using TranslationType = typename Superclass::TranslationType;
  using OffsetType = typename Superclass::OffsetType;
  using AngleType = typename Superclass::ScalarType;

  void SetParameters(const ParametersType &parameters) ITK_OVERRIDE;
  const ParametersType &GetParameters() const ITK_OVERRIDE;

  // const JacobianType & GetJacobian(const InputPointType  &point ) const;

  virtual void SetIdentity() ITK_OVERRIDE;

  // St a rotation center for the symmetry plane (easier for optimization)
  void SetRotationCenter(OutputPointType &rotCenter);

protected:
  SymmetryPlaneTransform();
  SymmetryPlaneTransform(const MatrixType &matrix,
                         const OutputPointType &offset);
  SymmetryPlaneTransform(unsigned int outputSpaceDims,
                         unsigned int paramsSpaceDims);

  virtual ~SymmetryPlaneTransform() {}

  void PrintSelf(std::ostream &os, itk::Indent indent) const ITK_OVERRIDE;

  void ComputeMatrix() ITK_OVERRIDE;

private:
  SymmetryPlaneTransform(const Self &); // purposely not implemented
  void operator=(const Self &);         // purposely not implemented

  // A symmetry plane is characterized by its spherical coordinates: theta, phi,
  // r
  ScalarType m_Theta;
  ScalarType m_Phi;
  ScalarType m_Distance;

  HomogeneousMatrixType m_CenterTransform, m_CenterTransformInv;
}; // class SymmetryPlaneTransform

} // end of namespace anima

#include "animaSymmetryPlaneTransform.hxx"
