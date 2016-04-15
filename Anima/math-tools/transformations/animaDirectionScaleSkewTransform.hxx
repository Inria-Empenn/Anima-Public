#pragma once
#include "animaDirectionScaleSkewTransform.h"

namespace anima
{

// Set Parameters
template <class TScalarType>
void
DirectionScaleSkewTransform<TScalarType>
::SetParameters(const ParametersType & parameters)
{
    itkDebugMacro( << "Setting parameters " << parameters );

    // Set angles with parameters
    m_LogScale = parameters[0];
    m_LogFirstSkew = parameters[1];
    m_LogSecondSkew = parameters[2];
    m_UniqueTranslation = parameters[3];

    this->ComputeMatrix();
    // Modified is always called since we just have a pointer to the
    // parameters and cannot know if the parameters have changed.
    this->Modified();

    itkDebugMacro(<<"After setting parameters ");
}


// Get Parameters
template <class TScalarType>
const typename DirectionScaleSkewTransform<TScalarType>::ParametersType &
DirectionScaleSkewTransform<TScalarType>
::GetParameters( void ) const
{
    this->m_Parameters[0] = m_LogScale;
    this->m_Parameters[1] = m_LogFirstSkew;
    this->m_Parameters[2] = m_LogSecondSkew;
    this->m_Parameters[3] = m_UniqueTranslation;

    return this->m_Parameters;
}

template <class TScalarType>
void
DirectionScaleSkewTransform<TScalarType>
::SetIdentity()
{
    Superclass::SetIdentity();
    m_LogScale = 0;
    m_UniqueTranslation = 0;
    m_LogFirstSkew = 0;
    m_LogSecondSkew = 0;

    this->ComputeMatrix();
}

template <class TScalarType>
void
DirectionScaleSkewTransform<TScalarType>
::SetGeometry(HomogeneousMatrixType &matrix)
{
    m_Geometry = matrix;
    m_GeometryInv = m_Geometry.GetInverse();
    this->ComputeMatrix();
}

template <class TScalarType>
void
DirectionScaleSkewTransform<TScalarType>
::SetLogScale(ScalarType scale, bool update)
{
    m_LogScale = scale;

    if (update)
        this->ComputeMatrix();
}

template <class TScalarType>
void
DirectionScaleSkewTransform<TScalarType>
::SetLogFirstSkew(ScalarType skew, bool update)
{
    m_LogFirstSkew = skew;

    if (update)
        this->ComputeMatrix();
}

template <class TScalarType>
void
DirectionScaleSkewTransform<TScalarType>
::SetLogSecondSkew(ScalarType skew, bool update)
{
    m_LogSecondSkew = skew;

    if (update)
        this->ComputeMatrix();
}

template <class TScalarType>
void
DirectionScaleSkewTransform<TScalarType>
::SetUniqueTranslation(ScalarType translation, bool update)
{
    m_UniqueTranslation = translation;

    if (update)
        this->ComputeMatrix();
}

template <class TScalarType>
void
DirectionScaleSkewTransform<TScalarType>
::SetUniqueDirection(unsigned int direction)
{
    m_UniqueDirection = direction;
    this->ComputeMatrix();
}

// Compute angles from the rotation matrix
template <class TScalarType>
void
DirectionScaleSkewTransform<TScalarType>
::ComputeMatrix(void)
{
    this->ComputeLogRepresentations();

    vnl_matrix <TScalarType> transformMatrix = anima::GetExponential(m_LogTransform);

    MatrixType linearMatrix;
    OffsetType off;

    for (unsigned int i = 0;i < 3;++i)
    {
        off[i] = transformMatrix(i,3);
        for (unsigned int j = 0;j < 3;++j)
            linearMatrix(i,j) = transformMatrix(i,j);
    }

    this->SetMatrix(linearMatrix);
    this->SetOffset(off);
}

template <class TScalarType>
void
DirectionScaleSkewTransform<TScalarType>
::ComputeLogRepresentations()
{
    m_LogTransform.fill(itk::NumericTraits<TScalarType>::Zero);

    unsigned int skewLocation = 0;
    if (m_UniqueDirection == skewLocation)
        skewLocation++;

    m_LogTransform(m_UniqueDirection,skewLocation) = m_LogFirstSkew;
    skewLocation++;
    if (m_UniqueDirection == skewLocation)
        skewLocation++;

    m_LogTransform(m_UniqueDirection,skewLocation) = m_LogSecondSkew;

    m_LogTransform(m_UniqueDirection,m_UniqueDirection) = m_LogScale;
    m_LogTransform(m_UniqueDirection,3) = m_UniqueTranslation;

    m_LogTransform = m_Geometry.GetVnlMatrix() * m_LogTransform * m_GeometryInv.GetVnlMatrix();

    m_LogVector.Fill(itk::NumericTraits<TScalarType>::Zero);

    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j <= 3;++j)
            m_LogVector[i * 4 + j] = m_LogTransform(i,j);
}

// Print self
template<class TScalarType>
void
DirectionScaleSkewTransform<TScalarType>::
PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os,indent);

    os << indent << "Direction scale transform parameters: Log scale=" << m_LogScale
       << " Log first skew=" << m_LogFirstSkew
       << " Log second skew=" << m_LogSecondSkew
       << " Log translation=" << m_UniqueTranslation
       << " Transform direction=" << m_UniqueDirection
       << " Geometry=" << m_Geometry
       << std::endl;
}

// Set Parameters
template <class TScalarType>
void
DirectionScaleTransform<TScalarType>
::SetParameters(const ParametersType & parameters)
{
    itkDebugMacro( << "Setting parameters " << parameters );

    // Set angles with parameters
    this->SetLogScale(parameters[0],false);
    this->SetLogFirstSkew(0,false);
    this->SetLogSecondSkew(0,false);
    this->SetUniqueTranslation(parameters[1],false);

    this->ComputeMatrix();
    // Modified is always called since we just have a pointer to the
    // parameters and cannot know if the parameters have changed.
    this->Modified();

    itkDebugMacro(<<"After setting parameters ");
}


// Get Parameters
template <class TScalarType>
const typename DirectionScaleTransform<TScalarType>::ParametersType &
DirectionScaleTransform<TScalarType>
::GetParameters( void ) const
{
    this->m_Parameters[0] = this->GetLogScale();
    this->m_Parameters[1] = this->GetUniqueTranslation();

    return this->m_Parameters;
}

// Set Parameters
template <class TScalarType>
void
DirectionTransform<TScalarType>
::SetParameters(const ParametersType & parameters)
{
    itkDebugMacro( << "Setting parameters " << parameters );

    // Set angles with parameters
    this->SetLogScale(0,false);
    this->SetLogFirstSkew(0,false);
    this->SetLogSecondSkew(0,false);
    this->SetUniqueTranslation(parameters[0],false);

    this->ComputeMatrix();
    // Modified is always called since we just have a pointer to the
    // parameters and cannot know if the parameters have changed.
    this->Modified();

    itkDebugMacro(<<"After setting parameters ");
}


// Get Parameters
template <class TScalarType>
const typename DirectionTransform<TScalarType>::ParametersType &
DirectionTransform<TScalarType>
::GetParameters() const
{
    this->m_Parameters[0] = this->GetUniqueTranslation();
    return this->m_Parameters;
}

} // end namespace anima
