#pragma once

#include <iostream>
#include <itkMatrixOffsetTransformBase.h>

namespace anima {

template <class TScalarType = double>
class LogRigid3DTransform
    : public itk::MatrixOffsetTransformBase<TScalarType, 3> {
public:
  /** Standard class typedefs. */
  using Self = LogRigid3DTransform;
  using Superclass = itk::MatrixOffsetTransformBase<TScalarType, 3>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** New macro for creation of through a Smart Pointer. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(LogRigid3DTransform, MatrixOffsetTransformBase);

  /** Dimension of the space. */
  itkStaticConstMacro(SpaceDimension, unsigned int, 3);
  itkStaticConstMacro(InputSpaceDimension, unsigned int, 3);
  itkStaticConstMacro(OutputSpaceDimension, unsigned int, 3);
  itkStaticConstMacro(ParametersDimension, unsigned int, 6);

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

  /** Set the transformation from a container of parameters
   * This is typically used by optimizers.
   * There are 6 parameters:
   *   0-3  angles
   *   3-5  translation
   **  */
  virtual void SetParameters(const ParametersType &parameters) ITK_OVERRIDE;
  virtual const ParametersType &GetParameters() const ITK_OVERRIDE;

  virtual const itk::Vector<TScalarType, 6> &GetLogVector() const {
    return m_LogVector;
  }
  virtual const vnl_matrix<TScalarType> &GetLogTransform() const {
    return m_LogTransform;
  }

  using AngleVectorType = itk::Vector<TScalarType, 3>;

  void SetAngles(const AngleVectorType &logAngle);
  itkGetConstReferenceMacro(Angles, AngleVectorType);

  virtual void SetIdentity() ITK_OVERRIDE;

protected:
  LogRigid3DTransform();
  virtual ~LogRigid3DTransform() {}

  void PrintSelf(std::ostream &os, itk::Indent indent) const ITK_OVERRIDE;

  void ComputeMatrix() ITK_OVERRIDE;
  void ComputeLogRepresentations();

private:
  LogRigid3DTransform(const Self &); // purposely not implemented
  void operator=(const Self &);      // purposely not implemented

  /**  Vector containing the angles. */
  AngleVectorType m_Angles;

  itk::Vector<TScalarType, 6> m_LogVector;
  vnl_matrix<TScalarType> m_LogTransform;

}; // class LogRigid3DTransform

} // end of namespace anima

#include "animaLogRigid3DTransform.hxx"
