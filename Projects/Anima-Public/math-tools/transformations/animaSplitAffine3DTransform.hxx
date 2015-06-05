#pragma once

#include "animaSplitAffine3DTransform.h"
#include <itkMath.h>

namespace anima
{

// Constructor with default arguments
template <class TScalarType>
SplitAffine3DTransform<TScalarType>
::SplitAffine3DTransform()
: Superclass(ParametersDimension)
{
    m_Angle.Fill( itk::NumericTraits<TScalarType>::Zero );
    m_LogScale.Fill( itk::NumericTraits<TScalarType>::One );
    m_Skew.Fill( itk::NumericTraits<TScalarType>::Zero );
}

// Set Parameters
template <class TScalarType>
void
SplitAffine3DTransform<TScalarType>
::SetParameters( const ParametersType & parameters)
{
    itkDebugMacro( << "Setting parameters " << parameters );

    // Transfer the versor part

    m_Angle[0] = parameters[0];
    m_Angle[1] = parameters[1];
    m_Angle[2] = parameters[2];

    // Matrix must be defined before translation so that offset can be computed
    // from translation
    m_LogScale[0] = parameters[6];
    m_LogScale[1] = parameters[7];
    m_LogScale[2] = parameters[8];

    m_Skew[0] = parameters[9];
    m_Skew[1] = parameters[10];
    m_Skew[2] = parameters[11];

    // Transfer the translation part
    TranslationType newTranslation;
    newTranslation[0] = parameters[3];
    newTranslation[1] = parameters[4];
    newTranslation[2] = parameters[5];

    this->SetTranslation(newTranslation);

    this->ComputeMatrix();
    this->ComputeOffset();

    // Modified is always called since we just have a pointer to the
    // parameters and cannot know if the parameters have changed.
    this->Modified();

    itkDebugMacro(<<"After setting parameters ");
}

//
// Get Parameters
//
// Parameters are ordered as:
//
// p[0:2] = euler angles (in radians)
// p[3:5] = translation components
// p[6:8] = Scale
// p[9:14] = Skew {xy, xz, yx, yz, zx, zy}
//

template <class TScalarType>
const typename SplitAffine3DTransform<TScalarType>::ParametersType &
SplitAffine3DTransform<TScalarType>
::GetParameters( void ) const
{
    itkDebugMacro( << "Getting parameters ");

    this->m_Parameters[0] = this->GetAngle()[0];
    this->m_Parameters[1] = this->GetAngle()[1];
    this->m_Parameters[2] = this->GetAngle()[2];

    this->m_Parameters[3] = this->GetTranslation()[0];
    this->m_Parameters[4] = this->GetTranslation()[1];
    this->m_Parameters[5] = this->GetTranslation()[2];

    this->m_Parameters[6] = this->GetLogScale()[0];
    this->m_Parameters[7] = this->GetLogScale()[1];
    this->m_Parameters[8] = this->GetLogScale()[2];

    this->m_Parameters[9] = this->GetSkew()[0];
    this->m_Parameters[10] = this->GetSkew()[1];
    this->m_Parameters[11] = this->GetSkew()[2];

    itkDebugMacro(<<"After getting parameters " << this->m_Parameters );

    return this->m_Parameters;
}

template <class TScalarType>
void
SplitAffine3DTransform<TScalarType>
::SetIdentity()
{
    m_Angle.Fill( itk::NumericTraits< ScaleVectorValueType >::Zero );
    m_LogScale.Fill( itk::NumericTraits< ScaleVectorValueType >::Zero );
    m_Skew.Fill( itk::NumericTraits< SkewVectorValueType >::Zero );

    TranslationType newTranslation;
    newTranslation.Fill( 0 );

    this->SetTranslation(newTranslation);

    this->ComputeMatrix();
    this->ComputeOffset();
}

template <class TScalarType>
void
SplitAffine3DTransform<TScalarType>
::SetAngle( const AngleVectorType & angle )
{
    m_Angle = angle;
    this->ComputeMatrix();
    this->ComputeOffset();
}

template <class TScalarType>
void
SplitAffine3DTransform<TScalarType>
::SetLogScale( const ScaleVectorType & scale )
{
    m_LogScale = scale;
    this->ComputeMatrix();
    this->ComputeOffset();
}

template <class TScalarType>
void
SplitAffine3DTransform<TScalarType>
::SetSkew( const SkewVectorType & skew )
{
    m_Skew = skew;
    this->ComputeMatrix();
    this->ComputeOffset();
}

// Compute the matrix
template <class TScalarType>
void
SplitAffine3DTransform<TScalarType>
::ComputeMatrix( void )
{
    MatrixType rotationMatrix = this->ComputeRotationMatrix();

    MatrixType scaleMatrix;
    scaleMatrix.SetIdentity();
    for (unsigned int i = 0;i < InputSpaceDimension;++i)
        scaleMatrix[i][i] = std::exp(m_LogScale[i]);

    MatrixType skewMatrix = this->ComputeSkewMatrix();

    this->SetVarMatrix ( rotationMatrix * scaleMatrix * skewMatrix );
}

template <class TScalarType>
typename SplitAffine3DTransform<TScalarType>::MatrixType
SplitAffine3DTransform<TScalarType>
::ComputeSkewMatrix()
{
    MatrixType skewMatrix;
    skewMatrix.SetIdentity();

    skewMatrix[0][1] = tan(m_Skew[0]);
    skewMatrix[0][2] = tan(m_Skew[1]);
    skewMatrix[1][2] = tan(m_Skew[2]);

    return skewMatrix;
}

template <class TScalarType>
typename SplitAffine3DTransform<TScalarType>::MatrixType
SplitAffine3DTransform<TScalarType>
::ComputeRotationMatrix()
{
    const ScalarType cx = vcl_cos(m_Angle[0]);
    const ScalarType sx = vcl_sin(m_Angle[0]);
    const ScalarType cy = vcl_cos(m_Angle[1]);
    const ScalarType sy = vcl_sin(m_Angle[1]);
    const ScalarType cz = vcl_cos(m_Angle[2]);
    const ScalarType sz = vcl_sin(m_Angle[2]);
    const ScalarType one = itk::NumericTraits< ScalarType >::One;
    const ScalarType zero = itk::NumericTraits< ScalarType >::Zero;

    MatrixType RotationX;
    RotationX[0][0]=one;  RotationX[0][1]=zero; RotationX[0][2]=zero;
    RotationX[1][0]=zero; RotationX[1][1]=cx;   RotationX[1][2]=-sx;
    RotationX[2][0]=zero; RotationX[2][1]=sx;   RotationX[2][2]=cx;


    MatrixType RotationY;
    RotationY[0][0]=cy;   RotationY[0][1]=zero; RotationY[0][2]=sy;
    RotationY[1][0]=zero; RotationY[1][1]=one;  RotationY[1][2]=zero;
    RotationY[2][0]=-sy;  RotationY[2][1]=zero; RotationY[2][2]=cy;


    MatrixType RotationZ;
    RotationZ[0][0]=cz;   RotationZ[0][1]=-sz;  RotationZ[0][2]=zero;
    RotationZ[1][0]=sz;   RotationZ[1][1]=cz;   RotationZ[1][2]=zero;
    RotationZ[2][0]=zero; RotationZ[2][1]=zero; RotationZ[2][2]=one;

    return RotationZ * RotationX * RotationY;
}

template <class TScalarType>
void
SplitAffine3DTransform<TScalarType>
::ComputeMatrixParameters( void )
{
    itkExceptionMacro( << "Setting the matrix of a SplitAffine3D transform is not supported at this time." );
}


// Print self
template<class TScalarType>
void
SplitAffine3DTransform<TScalarType>::
PrintSelf(std::ostream &os, itk::Indent indent) const
{

    Superclass::PrintSelf(os,indent);

    os << indent << "Angle:       " << m_Angle        << std::endl;
    os << indent << "Scale:       " << m_LogScale        << std::endl;
    os << indent << "Skew:        " << m_Skew         << std::endl;
}

} // end of namespace anima
