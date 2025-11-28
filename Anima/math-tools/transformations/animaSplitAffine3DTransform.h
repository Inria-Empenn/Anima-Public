#pragma once

#include <iostream>
#include <itkMatrixOffsetTransformBase.h>

namespace anima {

template <class TScalarType = double> // Data type for scalars:double or double
class SplitAffine3DTransform
    : public itk::MatrixOffsetTransformBase<TScalarType, 3, 3> {
public:
  /** Standard class typedefs. */
  using Self = SplitAffine3DTransform;
  using Superclass = itk::MatrixOffsetTransformBase<TScalarType, 3, 3>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** New macro for creation of through a Smart Pointer. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SplitAffine3DTransform, MatrixOffsetTransformBase);

  /** Dimension of parameters. */
  itkStaticConstMacro(InputSpaceDimension, unsigned int, 3);
  itkStaticConstMacro(OutputSpaceDimension, unsigned int, 3);
  itkStaticConstMacro(ParametersDimension, unsigned int, 12);

  /** Parameters Type   */
  using ParametersType = typename Superclass::ParametersType;
  using ScalarType = typename Superclass::ScalarType;
  using MatrixType = typename Superclass::MatrixType;
  using OffsetType = typename Superclass::OffsetType;
  using TranslationType = typename Superclass::TranslationType;

  /** Scale & Skew Vector Type. */
  using VectorType = itk::Vector<TScalarType, 3>;

  /** Set the transformation from a container of parameters
   * This is typically used by optimizers.
   * There are 12 parameters:
   *   0-2   euler angles (first matrix)
   *   3-5   translation
   *   6-8   Scale
   *   9-11  euler angles (second matrix)
   **  */
  virtual void SetParameters(const ParametersType &parameters) ITK_OVERRIDE;
  virtual const ParametersType &GetParameters() const ITK_OVERRIDE;

  void SetFirstAngle(const VectorType &angle);
  itkGetConstReferenceMacro(FirstAngle, VectorType);

  void SetLogScale(const VectorType &scale);
  itkGetConstReferenceMacro(LogScale, VectorType);

  void SetSecondAngle(const VectorType &angle);
  itkGetConstReferenceMacro(SecondAngle, VectorType);

  void SetIdentity() ITK_OVERRIDE;

protected:
  SplitAffine3DTransform();
  virtual ~SplitAffine3DTransform() {}

  void PrintSelf(std::ostream &os, itk::Indent indent) const ITK_OVERRIDE;

  MatrixType ComputeRotationMatrix();
  MatrixType ComputeSkewMatrix();

  /** Compute the components of the rotation matrix in the superclass. */
  void ComputeMatrix() ITK_OVERRIDE;
  void ComputeMatrixParameters() ITK_OVERRIDE;

private:
  SplitAffine3DTransform(const Self &); // purposely not implemented
  void operator=(const Self &);         // purposely not implemented

  /**  Vectors containing the euler angles. */
  VectorType m_FirstAngle, m_SecondAngle;

  /**  Vector containing the log-scale (exped when used). */
  VectorType m_LogScale;
};

} // end of namespace anima

#include "animaSplitAffine3DTransform.hxx"
