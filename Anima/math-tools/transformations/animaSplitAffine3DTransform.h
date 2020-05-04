#pragma once

#include <iostream>
#include <itkMatrixOffsetTransformBase.h>

namespace anima
{

template < class TScalarType=double >  // Data type for scalars:double or double
class SplitAffine3DTransform :
public itk::MatrixOffsetTransformBase< TScalarType, 3, 3 >
{
public:
    /** Standard class typedefs. */
    typedef SplitAffine3DTransform Self;
    typedef itk::MatrixOffsetTransformBase <TScalarType, 3, 3> Superclass;
    typedef itk::SmartPointer <Self> Pointer;
    typedef itk::SmartPointer <const Self> ConstPointer;

    /** New macro for creation of through a Smart Pointer. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(SplitAffine3DTransform, MatrixOffsetTransformBase)

    /** Dimension of parameters. */
    itkStaticConstMacro(InputSpaceDimension, unsigned int, 3);
    itkStaticConstMacro(OutputSpaceDimension, unsigned int, 3);
    itkStaticConstMacro(ParametersDimension, unsigned int, 12);

    /** Parameters Type   */
    typedef typename Superclass::ParametersType ParametersType;
    typedef typename Superclass::ScalarType ScalarType;
    typedef typename Superclass::MatrixType MatrixType;
    typedef typename Superclass::OffsetType OffsetType;
    typedef typename Superclass::TranslationType TranslationType;

    /** Scale & Skew Vector Type. */
    typedef itk::Vector<TScalarType, 3> VectorType;

    /** Set the transformation from a container of parameters
     * This is typically used by optimizers.
     * There are 12 parameters:
     *   0-2   euler angles (first matrix)
     *   3-5   translation
     *   6-8   Scale
     *   9-11  euler angles (second matrix)
     **  */
    virtual void SetParameters(const ParametersType & parameters) ITK_OVERRIDE;
    virtual const ParametersType& GetParameters() const ITK_OVERRIDE;

    void SetFirstAngle(const VectorType &angle);
    itkGetConstReferenceMacro(FirstAngle, VectorType)

    void SetLogScale(const VectorType &scale);
    itkGetConstReferenceMacro(LogScale, VectorType)

    void SetSecondAngle(const VectorType &angle);
    itkGetConstReferenceMacro(SecondAngle, VectorType)

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
    SplitAffine3DTransform(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    /**  Vectors containing the euler angles. */
    VectorType m_FirstAngle, m_SecondAngle;

    /**  Vector containing the log-scale (exped when used). */
    VectorType m_LogScale;
};

}  // end of namespace anima

#include "animaSplitAffine3DTransform.hxx"
