#pragma once

#include <iostream>
#include <itkMatrixOffsetTransformBase.h>

namespace anima
{

template < class TScalarType=double >    // Data type for scalars (float or double)
class SymmetryConstrainedTransform :
public itk::MatrixOffsetTransformBase< TScalarType, 3 >
{
public:
    /** Standard class typedefs. */
    typedef SymmetryConstrainedTransform                  Self;
    typedef itk::MatrixOffsetTransformBase< TScalarType, 3 >   Superclass;
    typedef itk::SmartPointer<Self>                Pointer;
    typedef itk::SmartPointer<const Self>          ConstPointer;

    /** New macro for creation of through a Smart Pointer. */
    itkNewMacro( Self );

    /** Run-time type information (and related methods). */
    itkTypeMacro( SymmetryConstrainedTransform, MatrixOffsetTransformBase );

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

    void SetParameters( const ParametersType & parameters );
    const ParametersType& GetParameters(void) const;

    //const JacobianType & GetJacobian(const InputPointType  &point ) const;

    virtual void SetIdentity(void);

    void SetReferenceSymmetryPlanes(Superclass *refSymmetryPlane, Superclass *floSymmetryPlane);

    void SetRotationCenter(OutputPointType &rotCenter);

protected:
    SymmetryConstrainedTransform();
    SymmetryConstrainedTransform(unsigned int outputSpaceDims,
                                 unsigned int paramsSpaceDims);

    virtual ~SymmetryConstrainedTransform(){}

    void PrintSelf(std::ostream &os, itk::Indent indent) const;

    void ComputeMatrix(void);

private:
    SymmetryConstrainedTransform(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    HomogeneousMatrixType m_RefSymmetryMatrix, m_FloSymmetryMatrix;

    ScalarType m_Angle;
    ScalarType m_TranslationY;
    ScalarType m_TranslationZ;

    HomogeneousMatrixType m_CenterTransform, m_CenterTransformInv;
}; //class SymmetryConstrainedTransform


}  // end of namespace anima

#include "animaSymmetryConstrainedTransform.hxx"
