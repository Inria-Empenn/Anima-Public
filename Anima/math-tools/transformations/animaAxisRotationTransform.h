#pragma once

#include <iostream>
#include <itkMatrixOffsetTransformBase.h>

namespace anima {

template <class TScalarType =
              double> // Data type for scalars (double or double)
class AxisRotationTransform
    : public itk::MatrixOffsetTransformBase<TScalarType, 3> {
public:
  /** Standard class typedefs. */
  using Self = AxisRotationTransform;
  using Superclass = itk::MatrixOffsetTransformBase<TScalarType, 3>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** New macro for creation of through a Smart Pointer. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(AxisRotationTransform, MatrixOffsetTransformBase);

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

  virtual void SetIdentity() ITK_OVERRIDE;

  void SetReferencePlaneNormal(InputVectorType &refPlaneNormal);

protected:
  AxisRotationTransform();
  AxisRotationTransform(unsigned int outputSpaceDims,
                        unsigned int paramsSpaceDims);

  virtual ~AxisRotationTransform() {}

  void PrintSelf(std::ostream &os, itk::Indent indent) const ITK_OVERRIDE;

  void ComputeMatrix() ITK_OVERRIDE;

private:
  AxisRotationTransform(const Self &); // purposely not implemented
  void operator=(const Self &);        // purposely not implemented

  HomogeneousMatrixType m_ReferencePlane, m_InverseReferencePlane;

  ScalarType m_Angle;
  ScalarType m_TranslationY;
  ScalarType m_TranslationZ;
}; // class AxisRotationTransform

} // end of namespace anima

#include "animaAxisRotationTransform.hxx"
