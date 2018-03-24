#pragma once

namespace anima
{

// Constructor with default arguments
template <class TScalarType>
SymmetryPlaneTransform<TScalarType>
::SymmetryPlaneTransform():
Superclass(ParametersDimension)
{
    m_Theta = m_Phi = m_Distance = itk::NumericTraits<ScalarType>::Zero;

    m_CenterTransform.SetIdentity();
    m_CenterTransformInv.SetIdentity();

    this->ComputeMatrix();
}

// Constructor with default arguments
template <class TScalarType>
SymmetryPlaneTransform<TScalarType>
::SymmetryPlaneTransform(const MatrixType & matrix,
                         const OutputPointType & offset)
{
    this->SetMatrix(matrix);

    OffsetType off;
    off[0] = offset[0];
    off[1] = offset[1];
    off[2] = offset[2];
    this->SetOffset(off);
}


// Constructor with arguments
template <class TScalarType>
SymmetryPlaneTransform<TScalarType>
::SymmetryPlaneTransform(unsigned int spaceDimension,
                         unsigned int parametersDimension):
Superclass(spaceDimension, parametersDimension)
{
    m_Theta = m_Phi = m_Distance = itk::NumericTraits<ScalarType>::Zero;

    m_CenterTransform.SetIdentity();
    m_CenterTransformInv.SetIdentity();

    this->ComputeMatrix();
}

template <class TScalarType>
void
SymmetryPlaneTransform<TScalarType>
::SetRotationCenter(OutputPointType &rotCenter)
{
    m_CenterTransform.SetIdentity();
    m_CenterTransformInv.SetIdentity();

    for (unsigned int i = 0;i < 3;++i)
    {
        m_CenterTransform(i,3) = rotCenter[i];
        m_CenterTransformInv(i,3) = - rotCenter[i];
    }

    this->ComputeMatrix();
    // Modified is always called since we just have a pointer to the
    // parameters and cannot know if the parameters have changed.
    this->Modified();
}

// Set Parameters
template <class TScalarType>
void
SymmetryPlaneTransform<TScalarType>
::SetParameters( const ParametersType & parameters )
{
    itkDebugMacro( << "Setting parameters " << parameters );

    // Set angles with parameters
    m_Theta = parameters[0];
    m_Phi = parameters[1];
    m_Distance = parameters[2];

    this->ComputeMatrix();
    // Modified is always called since we just have a pointer to the
    // parameters and cannot know if the parameters have changed.
    this->Modified();

    itkDebugMacro(<<"After setting parameters ");
}


// Get Parameters
template <class TScalarType>
const typename SymmetryPlaneTransform<TScalarType>::ParametersType &
SymmetryPlaneTransform<TScalarType>
::GetParameters( void ) const
{
    this->m_Parameters[0] = m_Theta;
    this->m_Parameters[1] = m_Phi;
    this->m_Parameters[2] = m_Distance;

    return this->m_Parameters;
}

// Compose
template <class TScalarType>
void
SymmetryPlaneTransform<TScalarType>
::SetIdentity(void)
{
    Superclass::SetIdentity();
    m_Theta = 0;
    m_Phi = 0;
    m_Distance = 0;

    this->ComputeMatrix();
}

// Compute angles from the rotation matrix
template <class TScalarType>
void
SymmetryPlaneTransform<TScalarType>
::ComputeMatrix(void)
{
    double cosa, sina, cosb, sinb;

    cosa = std::cos(m_Theta);
    sina = std::sin(m_Theta);
    cosb = std::cos(m_Phi);
    sinb = std::sin(m_Phi);

    double d = m_Distance;

    HomogeneousMatrixType rotationMatrix;
    rotationMatrix.SetIdentity();

    rotationMatrix[0][0] = 1-2*cosb*cosb*cosa*cosa;
    rotationMatrix[0][1] = -2*cosb*sinb*cosa*cosa;
    rotationMatrix[0][2] = -2*cosb*cosa*sina;

    rotationMatrix[1][0] = rotationMatrix[0][1];
    rotationMatrix[1][1] = 1-2*sinb*sinb*cosa*cosa;
    rotationMatrix[1][2] = -2*sinb*cosa*sina;

    rotationMatrix[2][0] = rotationMatrix[0][2];
    rotationMatrix[2][1] = rotationMatrix[1][2];
    rotationMatrix[2][2] = 1-2*sina*sina;

    rotationMatrix[0][3] = 2*d*cosb*cosa;
    rotationMatrix[1][3] = 2*d*sinb*cosa;
    rotationMatrix[2][3] = 2*d*sina;

    rotationMatrix = m_CenterTransform * rotationMatrix * m_CenterTransformInv;

    MatrixType linearMatrix;
    OffsetType off;

    for (unsigned int i = 0;i < 3;++i)
    {
        off[i] = rotationMatrix(i,3);
        for (unsigned int j = 0;j < 3;++j)
            linearMatrix(i,j) = rotationMatrix(i,j);
    }

    this->SetMatrix(linearMatrix);
    this->SetOffset(off);
}

// Print self
template<class TScalarType>
void
SymmetryPlaneTransform<TScalarType>::
PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os,indent);

    os << indent << "Symmetry plane parameters: Theta=" << m_Theta
    << " Phi=" << m_Phi
    << " Distance=" << m_Distance
    << " Center transform=" << m_CenterTransform
    << std::endl;
}

} // end of namespace anima
