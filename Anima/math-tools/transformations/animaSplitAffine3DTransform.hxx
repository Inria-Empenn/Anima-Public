#pragma once

#include "animaSplitAffine3DTransform.h"
#include <animaMatrixLogExp.h>

namespace anima
{

// Constructor with default arguments
template <class TScalarType>
SplitAffine3DTransform<TScalarType>
::SplitAffine3DTransform()
: Superclass(ParametersDimension)
{
    m_FirstAngle.Fill(itk::NumericTraits<TScalarType>::Zero);
    m_SecondAngle.Fill(itk::NumericTraits<TScalarType>::Zero);
    m_LogScale.Fill(itk::NumericTraits<TScalarType>::Zero);
}

// Set Parameters
template <class TScalarType>
void
SplitAffine3DTransform<TScalarType>
::SetParameters( const ParametersType & parameters)
{
    m_FirstAngle[0] = parameters[0];
    m_FirstAngle[1] = parameters[1];
    m_FirstAngle[2] = parameters[2];

    // Matrix must be defined before translation so that offset can be computed
    // from translation
    m_LogScale[0] = parameters[6];
    m_LogScale[1] = parameters[7];
    m_LogScale[2] = parameters[8];

    m_SecondAngle[0] = parameters[9];
    m_SecondAngle[1] = parameters[10];
    m_SecondAngle[2] = parameters[11];

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
}

//
// Get Parameters
//
template <class TScalarType>
const typename SplitAffine3DTransform<TScalarType>::ParametersType &
SplitAffine3DTransform<TScalarType>
::GetParameters( void ) const
{
    this->m_Parameters[0] = this->GetFirstAngle()[0];
    this->m_Parameters[1] = this->GetFirstAngle()[1];
    this->m_Parameters[2] = this->GetFirstAngle()[2];

    this->m_Parameters[3] = this->GetTranslation()[0];
    this->m_Parameters[4] = this->GetTranslation()[1];
    this->m_Parameters[5] = this->GetTranslation()[2];

    this->m_Parameters[6] = this->GetLogScale()[0];
    this->m_Parameters[7] = this->GetLogScale()[1];
    this->m_Parameters[8] = this->GetLogScale()[2];

    this->m_Parameters[9] = this->GetSecondAngle()[0];
    this->m_Parameters[10] = this->GetSecondAngle()[1];
    this->m_Parameters[11] = this->GetSecondAngle()[2];

    return this->m_Parameters;
}

template <class TScalarType>
void
SplitAffine3DTransform<TScalarType>
::SetIdentity()
{
    m_FirstAngle.Fill(itk::NumericTraits <TScalarType>::Zero);
    m_LogScale.Fill(itk::NumericTraits <TScalarType>::Zero);
    m_SecondAngle.Fill(itk::NumericTraits <TScalarType>::Zero);

    TranslationType newTranslation;
    newTranslation.Fill(0.0);

    this->SetTranslation(newTranslation);

    this->ComputeMatrix();
    this->ComputeOffset();
}

template <class TScalarType>
void
SplitAffine3DTransform<TScalarType>
::SetFirstAngle(const VectorType &angle)
{
    m_FirstAngle = angle;
    this->ComputeMatrix();
    this->ComputeOffset();
}

template <class TScalarType>
void
SplitAffine3DTransform<TScalarType>
::SetLogScale(const VectorType &scale)
{
    m_LogScale = scale;
    this->ComputeMatrix();
    this->ComputeOffset();
}

template <class TScalarType>
void
SplitAffine3DTransform<TScalarType>
::SetSecondAngle(const VectorType &angle)
{
    m_SecondAngle = angle;
    this->ComputeMatrix();
    this->ComputeOffset();
}

// Compute the matrix
template <class TScalarType>
void
SplitAffine3DTransform<TScalarType>
::ComputeMatrix()
{
    std::vector <TScalarType> angles(3, 0.0);
    for (unsigned int i = 0;i < InputSpaceDimension;++i)
        angles[i] = m_FirstAngle[i];

    vnl_matrix <TScalarType> firstRotationMatrix;
    anima::Get3DRotationExponential(angles,firstRotationMatrix);

    vnl_matrix <TScalarType> scaleMatrix;
    scaleMatrix.set_identity();
    for (unsigned int i = 0;i < InputSpaceDimension;++i)
        scaleMatrix(i,i) = std::exp(m_LogScale[i]);

    for (unsigned int i = 0;i < InputSpaceDimension;++i)
        angles[i] = m_SecondAngle[i];

    vnl_matrix <TScalarType> secondRotationMatrix;
    anima::Get3DRotationExponential(angles,secondRotationMatrix);

    MatrixType varMatrix = firstRotationMatrix * scaleMatrix * secondRotationMatrix;
    this->SetVarMatrix (varMatrix);
}

template <class TScalarType>
void
SplitAffine3DTransform<TScalarType>
::ComputeMatrixParameters()
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

    os << indent << "First angle: " << m_FirstAngle << std::endl;
    os << indent << "Scale: " << m_LogScale << std::endl;
    os << indent << "Second angle: " << m_SecondAngle << std::endl;
}

} // end of namespace anima
