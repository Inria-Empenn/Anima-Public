#pragma once

#include <iostream>
#include <itkMatrixOffsetTransformBase.h>

namespace anima
{

template < class TScalarType=double >
class LogRigid3DTransform :
public itk::MatrixOffsetTransformBase< TScalarType, 3 >
{
public:
    /** Standard class typedefs. */
    typedef LogRigid3DTransform                 Self;
    typedef itk::MatrixOffsetTransformBase <TScalarType, 3>   Superclass;
    typedef itk::SmartPointer<Self>                Pointer;
    typedef itk::SmartPointer<const Self>          ConstPointer;

    /** New macro for creation of through a Smart Pointer. */
    itkNewMacro( Self );

    /** Run-time type information (and related methods). */
    itkTypeMacro( LogRigid3DTransform, MatrixOffsetTransformBase );

    /** Dimension of the space. */
    itkStaticConstMacro(SpaceDimension, unsigned int, 3);
    itkStaticConstMacro(InputSpaceDimension, unsigned int, 3);
    itkStaticConstMacro(OutputSpaceDimension, unsigned int, 3);
    itkStaticConstMacro(ParametersDimension, unsigned int, 6);

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
    typedef itk::Matrix <ScalarType,4,4> HomogeneousMatrixType;
    typedef typename Superclass::InverseMatrixType          InverseMatrixType;
    typedef typename Superclass::CenterType                 CenterType;
    typedef typename Superclass::TranslationType            TranslationType;
    typedef typename Superclass::OffsetType                 OffsetType;

    /** Set the transformation from a container of parameters
     * This is typically used by optimizers.
     * There are 6 parameters:
     *   0-3  angles
     *   3-5  translation
     **  */
    virtual void SetParameters( const ParametersType & parameters );
    virtual const ParametersType& GetParameters() const;

    virtual const itk::Vector <TScalarType,6>& GetLogVector() const {return m_LogVector;}
    virtual const vnl_matrix <TScalarType>& GetLogTransform() const {return m_LogTransform;}

    typedef itk::Vector<TScalarType, 3> AngleVectorType;

    void SetAngles (const AngleVectorType & logAngle);
    itkGetConstReferenceMacro (Angles, AngleVectorType);

    virtual void SetIdentity();

protected:
    LogRigid3DTransform();
    virtual ~LogRigid3DTransform(){}

    void PrintSelf(std::ostream &os, itk::Indent indent) const;

    void ComputeMatrix();
    void ComputeLogRepresentations();

private:
    LogRigid3DTransform(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    /**  Vector containing the angles. */
    AngleVectorType m_Angles;

    itk::Vector <TScalarType,6> m_LogVector;
    vnl_matrix <TScalarType> m_LogTransform;

}; //class LogRigid3DTransform


}  // end of namespace anima

#include "animaLogRigid3DTransform.hxx"
