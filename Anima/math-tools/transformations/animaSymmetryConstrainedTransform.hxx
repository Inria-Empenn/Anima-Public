#pragma once

namespace anima
{

// Constructor with default arguments
template <class TScalarType>
SymmetryConstrainedTransform<TScalarType>
::SymmetryConstrainedTransform():
Superclass(ParametersDimension)
{
    m_Angle = m_TranslationY = m_TranslationZ = itk::NumericTraits<ScalarType>::Zero;
    m_RefSymmetryMatrix.SetIdentity();
    m_FloSymmetryMatrix.SetIdentity();

    m_CenterTransform.SetIdentity();
    m_CenterTransformInv.SetIdentity();

    this->ComputeMatrix();
}

// Constructor with arguments
template <class TScalarType>
SymmetryConstrainedTransform<TScalarType>
::SymmetryConstrainedTransform(unsigned int spaceDimension,
                         unsigned int parametersDimension):
Superclass(spaceDimension, parametersDimension)
{
    m_Angle = m_TranslationY = m_TranslationZ = itk::NumericTraits<ScalarType>::Zero;
    m_RefSymmetryMatrix.SetIdentity();
    m_FloSymmetryMatrix.SetIdentity();

    m_CenterTransform.SetIdentity();
    m_CenterTransformInv.SetIdentity();

    this->ComputeMatrix();
}

template <class TScalarType>
void
SymmetryConstrainedTransform<TScalarType>
::SetRotationCenter(OutputPointType &rotCenter)
{
    m_CenterTransform.SetIdentity();
    m_CenterTransformInv.SetIdentity();

    for (unsigned int i = 0;i < 3;++i)
    {
        m_CenterTransform(i,3) = rotCenter[i];
        m_CenterTransformInv(i,3) = - rotCenter[i];
    }
}

// Set Parameters
template <class TScalarType>
void
SymmetryConstrainedTransform<TScalarType>
::SetParameters( const ParametersType & parameters )
{
    itkDebugMacro( << "Setting parameters " << parameters );

    m_Angle = parameters[0];
    m_TranslationY = parameters[1];
    m_TranslationZ = parameters[2];

    this->ComputeMatrix();
    // Modified is always called since we just have a pointer to the
    // parameters and cannot know if the parameters have changed.
    this->Modified();

    itkDebugMacro(<<"After setting parameters ");
}


// Get Parameters
template <class TScalarType>
const typename SymmetryConstrainedTransform<TScalarType>::ParametersType &
SymmetryConstrainedTransform<TScalarType>
::GetParameters( void ) const
{
    this->m_Parameters[0] = m_Angle;
    this->m_Parameters[1] = m_TranslationY;
    this->m_Parameters[2] = m_TranslationZ;

    return this->m_Parameters;
}

// Compose
template <class TScalarType>
void
SymmetryConstrainedTransform<TScalarType>
::SetIdentity(void)
{
    Superclass::SetIdentity();
    m_Angle = m_TranslationY = m_TranslationZ = itk::NumericTraits<ScalarType>::Zero;

    this->ComputeMatrix();
}

// Compute angles from the rotation matrix
template <class TScalarType>
void
SymmetryConstrainedTransform<TScalarType>
::ComputeMatrix(void)
{
    double cosa, sina;

    cosa = cos(m_Angle);
    sina = sin(m_Angle);

    HomogeneousMatrixType insideMatrix, resultMatrix;
    insideMatrix.SetIdentity();

    insideMatrix[1][1] = cosa;
    insideMatrix[1][2] = -sina;

    insideMatrix[2][1] = sina;
    insideMatrix[2][2] = cosa;

    insideMatrix[1][3] = m_TranslationY;
    insideMatrix[2][3] = m_TranslationZ;

    insideMatrix = m_CenterTransform * insideMatrix * m_CenterTransformInv;
    resultMatrix = m_FloSymmetryMatrix * insideMatrix * m_RefSymmetryMatrix;

    MatrixType linearMatrix;
    OffsetType off;

    for (unsigned int i = 0;i < 3;++i)
    {
        off[i] = resultMatrix(i,3);
        for (unsigned int j = 0;j < 3;++j)
            linearMatrix(i,j) = resultMatrix(i,j);
    }

    this->SetMatrix(linearMatrix);
    this->SetOffset(off);
}

template<class TScalarType>
void
SymmetryConstrainedTransform<TScalarType>::
SetReferenceSymmetryPlanes(Superclass *refSymmetryPlane, Superclass *floSymmetryPlane)
{
    // The input transforms are matrix transforms realigning each image onto the mid-sagittal plane

    m_RefSymmetryMatrix.SetIdentity();

    MatrixType tmpMatrix;
    OffsetType tmpOffset;

    tmpMatrix = refSymmetryPlane->GetMatrix();
    tmpOffset = refSymmetryPlane->GetOffset();

    for (unsigned int i = 0;i < 3;++i)
    {
        m_RefSymmetryMatrix(i,3) = tmpOffset[i];
        for (unsigned int j = 0;j < 3;++j)
            m_RefSymmetryMatrix(i,j) = tmpMatrix(i,j);
    }

    m_RefSymmetryMatrix = m_RefSymmetryMatrix.GetInverse();

    tmpMatrix = floSymmetryPlane->GetMatrix();
    tmpOffset = floSymmetryPlane->GetOffset();

    for (unsigned int i = 0;i < 3;++i)
    {
        m_FloSymmetryMatrix(i,3) = tmpOffset[i];
        for (unsigned int j = 0;j < 3;++j)
            m_FloSymmetryMatrix(i,j) = tmpMatrix(i,j);
    }
}

// Print self
template<class TScalarType>
void
SymmetryConstrainedTransform<TScalarType>::
PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os,indent);

    os << indent << "Symmetry plane parameters: Angle=" << m_Angle
    << " TranslationY=" << m_TranslationY
    << " TranslationZ=" << m_TranslationZ
    << std::endl;
}

} // end of namespace anima
