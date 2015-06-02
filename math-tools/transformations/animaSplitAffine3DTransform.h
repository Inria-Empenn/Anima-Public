#pragma once

#include <iostream>
#include <itkMatrixOffsetTransformBase.h>

namespace anima
{

template < class TScalarType=double >  // Data type for scalars:float or double
class SplitAffine3DTransform :
public itk::MatrixOffsetTransformBase< TScalarType, 3, 3 >
{
public:
    /** Standard class typedefs. */
    typedef SplitAffine3DTransform              Self;
    typedef itk::MatrixOffsetTransformBase< TScalarType, 3, 3 >   Superclass;
    typedef itk::SmartPointer<Self>                      Pointer;
    typedef itk::SmartPointer<const Self>                ConstPointer;

    /** New macro for creation of through a Smart Pointer. */
    itkNewMacro( Self );

    /** Run-time type information (and related methods). */
    itkTypeMacro( SplitAffine3DTransform, MatrixOffsetTransformBase );

    /** Dimension of parameters. */
    itkStaticConstMacro(InputSpaceDimension, unsigned int, 3);
    itkStaticConstMacro(OutputSpaceDimension, unsigned int, 3);
    itkStaticConstMacro(ParametersDimension, unsigned int, 12);

    /** Parameters Type   */
    typedef typename Superclass::ParametersType         ParametersType;
    typedef typename Superclass::JacobianType           JacobianType;
    typedef typename Superclass::ScalarType             ScalarType;
    typedef typename Superclass::InputPointType         InputPointType;
    typedef typename Superclass::OutputPointType        OutputPointType;
    typedef typename Superclass::InputVectorType        InputVectorType;
    typedef typename Superclass::OutputVectorType       OutputVectorType;
    typedef typename Superclass::InputVnlVectorType     InputVnlVectorType;
    typedef typename Superclass::OutputVnlVectorType    OutputVnlVectorType;
    typedef typename Superclass::InputCovariantVectorType
    InputCovariantVectorType;
    typedef typename Superclass::OutputCovariantVectorType
    OutputCovariantVectorType;
    typedef typename Superclass::MatrixType             MatrixType;
    typedef typename Superclass::InverseMatrixType      InverseMatrixType;
    typedef typename Superclass::CenterType             CenterType;
    typedef typename Superclass::OffsetType             OffsetType;
    typedef typename Superclass::TranslationType        TranslationType;

    /** Scale & Skew Vector Type. */
    typedef itk::Vector<TScalarType, 3>                      AngleVectorType;
    typedef itk::Vector<TScalarType, 3>                      ScaleVectorType;
    typedef itk::Vector<TScalarType, 3>                     SkewVectorType;

    typedef typename AngleVectorType::ValueType         AngleVectorValueType;
    typedef typename ScaleVectorType::ValueType         ScaleVectorValueType;
    typedef typename SkewVectorType::ValueType          SkewVectorValueType;
    typedef typename TranslationType::ValueType         TranslationValueType;

    /** Set the transformation from a container of parameters
     * This is typically used by optimizers.
     * There are 12 parameters:
     *   0-2   euler angles
     *   3-5   translation
     *   6-8   Scale
     *   9-11  Skew
     **  */
    virtual void SetParameters( const ParametersType & parameters );
    virtual const ParametersType& GetParameters(void) const;

    void SetAngle( const AngleVectorType & angle );
    itkGetConstReferenceMacro( Angle, AngleVectorType );

    void SetLogScale( const ScaleVectorType & scale );
    itkGetConstReferenceMacro( LogScale, ScaleVectorType );

    void SetSkew( const SkewVectorType & skew );
    itkGetConstReferenceMacro( Skew, SkewVectorType );

    void SetIdentity();

protected:
    SplitAffine3DTransform();
    virtual ~SplitAffine3DTransform(){}

    void PrintSelf(std::ostream &os, itk::Indent indent) const;

    MatrixType ComputeRotationMatrix();
    MatrixType ComputeSkewMatrix();

    /** Compute the components of the rotation matrix in the superclass. */
    void ComputeMatrix(void);
    void ComputeMatrixParameters(void);

private:
    SplitAffine3DTransform(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    /**  Vector containing the euler angles. */
    AngleVectorType          m_Angle;

    /**  Vector containing the log-scale (exped when used). */
    ScaleVectorType          m_LogScale;

    /**  Vector containing the skews (skews are angles, take the tangent to compute the matrix) */
    SkewVectorType           m_Skew;

}; //class SplitAffine3DTransform


}  // end of namespace anima

#include "animaSplitAffine3DTransform.hxx"
