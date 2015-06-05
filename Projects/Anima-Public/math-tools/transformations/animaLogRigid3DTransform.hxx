#pragma once

#include <animaMatrixLogExp.h>

namespace anima
{

// Constructor with default arguments
template <class TScalarType>
LogRigid3DTransform<TScalarType>
::LogRigid3DTransform():
Superclass(ParametersDimension)
{
    m_Angles.Fill(itk::NumericTraits<TScalarType>::Zero);
    m_LogVector.Fill(itk::NumericTraits<TScalarType>::Zero);
    m_LogTransform.set_size(4,4);

    this->ComputeMatrix();
}

// Set Parameters
template <class TScalarType>
void
LogRigid3DTransform<TScalarType>
::SetParameters( const ParametersType & parameters )
{
    itkDebugMacro( << "Setting parameters " << parameters );

    m_Angles[0] = parameters[0];
    m_Angles[1] = parameters[1];
    m_Angles[2] = parameters[2];

    TranslationType newTranslation;
    newTranslation[0] = parameters[3];
    newTranslation[1] = parameters[4];
    newTranslation[2] = parameters[5];

    this->SetTranslation(newTranslation);

    this->ComputeMatrix();

    // Modified is always called since we just have a pointer to the
    // parameters and cannot know if the parameters have changed.
    this->Modified();

    itkDebugMacro(<<"After setting parameters ");
}


// Get Parameters
template <class TScalarType>
const typename LogRigid3DTransform<TScalarType>::ParametersType &
LogRigid3DTransform<TScalarType>
::GetParameters() const
{
    this->m_Parameters[0] = m_Angles[0];
    this->m_Parameters[1] = m_Angles[1];
    this->m_Parameters[2] = m_Angles[2];

    this->m_Parameters[3] = this->GetTranslation()[0];
    this->m_Parameters[4] = this->GetTranslation()[1];
    this->m_Parameters[5] = this->GetTranslation()[2];

    return this->m_Parameters;
}

// Compose
template <class TScalarType>
void
LogRigid3DTransform<TScalarType>
::SetIdentity()
{
    Superclass::SetIdentity();
    m_Angles.Fill(itk::NumericTraits<TScalarType>::Zero);

    this->ComputeMatrix();
}

// Compute angles from the rotation matrix
template <class TScalarType>
void
LogRigid3DTransform<TScalarType>
::ComputeMatrix()
{
    this->ComputeLogRepresentations();

    vnl_matrix <TScalarType> transformMatrix = GetExponential(m_LogTransform);

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
LogRigid3DTransform<TScalarType>
::ComputeLogRepresentations()
{
    m_LogVector.Fill(itk::NumericTraits<TScalarType>::Zero);

    for (unsigned int i = 0;i < 3;++i)
        m_LogVector[i] = m_Angles[i];

    m_LogTransform.fill(itk::NumericTraits<TScalarType>::Zero);

    unsigned int pos = 0;
    for (unsigned int i = 0;i < 3;++i)
    {
        for (unsigned int j = i+1;j < 3;++j)
        {
            m_LogTransform(i,j) = m_Angles[pos];
            m_LogTransform(j,i) = - m_Angles[pos];
            ++pos;
        }
    }

    for (unsigned int i = 0;i < 3;++i)
    {
        m_LogVector[i+3] = this->GetTranslation()[i];
        for (unsigned int j = 0;j < 3;++j)
            m_LogVector[i+3] -= m_LogTransform(i,j) * this->GetCenter()[j];
    }

    for (unsigned int i = 0;i < 3;++i)
        m_LogTransform(i,3) = m_LogVector[i+3];
}

// Print self
template<class TScalarType>
void
LogRigid3DTransform<TScalarType>::
PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os,indent);

    os << indent << "Log-transform parameters: Angles=" << m_Angles << std::endl;
}

} // end of namespace anima
