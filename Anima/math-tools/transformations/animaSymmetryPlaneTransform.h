#pragma once

#include <iostream>
#include <itkMatrixOffsetTransformBase.h>

namespace anima
{

template < class TScalarType=double >    // Data type for scalars (float or double)
class SymmetryPlaneTransform :
public itk::MatrixOffsetTransformBase< TScalarType, 3 >
{
public:
    /** Standard class typedefs. */
    typedef SymmetryPlaneTransform                  Self;
    typedef itk::MatrixOffsetTransformBase< TScalarType, 3 >   Superclass;
    typedef itk::SmartPointer<Self>                Pointer;
    typedef itk::SmartPointer<const Self>          ConstPointer;

    /** New macro for creation of through a Smart Pointer. */
    itkNewMacro( Self );

    /** Run-time type information (and related methods). */
    itkTypeMacro( SymmetryPlaneTransform, MatrixOffsetTransformBase );

    /** Dimension of the space. */
    itkStaticConstMacro(SpaceDimension, unsigned int, 3);
    itkStaticConstMacro(InputSpaceDimension, unsigned int, 3);
    itkStaticConstMacro(OutputSpaceDimension, unsigned int, 3);
    itkStaticConstMacro(ParametersDimension, unsigned int, 3);

    //Check which typedefs are useful
    typedef typename Superclass::ParametersType             ParametersType;
    typedef typename Superclass::ParametersValueType        ParametersValueType;
    typedef typename Superclass::JacobianType               JacobianType;
    typedef typename Superclass::ScalarType                 ScalarType;
    typedef typename Superclass::InputVectorType            InputVectorType;
    typedef typename Superclass::OutputVectorType           OutputVectorType;
    typedef typename Superclass::InputCovariantVectorType   InputCovariantVectorType;
    typedef typename Superclass::OutputCovariantVectorType  OutputCovariantVectorType;
    typedef typename Superclass::InputVnlVectorType         InputVnlVectorType;
    typedef typename Superclass::OutputVnlVectorType        OutputVnlVectorType;
    typedef typename Superclass::InputPointType             InputPointType;
    typedef typename Superclass::OutputPointType            OutputPointType;
    typedef typename Superclass::MatrixType                 MatrixType;
    typedef itk::Matrix <ScalarType,4,4>                    HomogeneousMatrixType;
    typedef typename Superclass::InverseMatrixType          InverseMatrixType;
    typedef typename Superclass::CenterType                 CenterType;
    typedef typename Superclass::TranslationType            TranslationType;
    typedef typename Superclass::OffsetType                 OffsetType;
    typedef typename Superclass::ScalarType                 AngleType;

    void SetParameters(const ParametersType & parameters) ITK_OVERRIDE;
    const ParametersType& GetParameters() const ITK_OVERRIDE;

    //const JacobianType & GetJacobian(const InputPointType  &point ) const;

    virtual void SetIdentity() ITK_OVERRIDE;

    // St a rotation center for the symmetry plane (easier for optimization)
    void SetRotationCenter(OutputPointType &rotCenter);

protected:
    SymmetryPlaneTransform();
    SymmetryPlaneTransform(const MatrixType & matrix,
                           const OutputPointType & offset);
    SymmetryPlaneTransform(unsigned int outputSpaceDims,
                           unsigned int paramsSpaceDims);

    virtual ~SymmetryPlaneTransform(){}

    void PrintSelf(std::ostream &os, itk::Indent indent) const ITK_OVERRIDE;

    void ComputeMatrix() ITK_OVERRIDE;

private:
    SymmetryPlaneTransform(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    // A symmetry plane is characterized by its spherical coordinates: theta, phi, r
    ScalarType  m_Theta;
    ScalarType  m_Phi;
    ScalarType  m_Distance;

    HomogeneousMatrixType m_CenterTransform, m_CenterTransformInv;
}; //class SymmetryPlaneTransform

}  // end of namespace anima

#include "animaSymmetryPlaneTransform.hxx"
