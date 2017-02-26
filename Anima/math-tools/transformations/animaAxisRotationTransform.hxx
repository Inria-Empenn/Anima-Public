#pragma once
#include <animaAxisRotationTransform.h>
#include <animaMatrixOperations.h>

namespace anima
{

// Constructor with default arguments
template <class TScalarType>
AxisRotationTransform<TScalarType>
::AxisRotationTransform():
    Superclass(ParametersDimension)
{
    m_Angle = m_TranslationY = m_TranslationZ = itk::NumericTraits<ScalarType>::Zero;
    m_ReferencePlane.SetIdentity();
    m_InverseReferencePlane.SetIdentity();

    this->ComputeMatrix();
}

// Constructor with arguments
template <class TScalarType>
AxisRotationTransform<TScalarType>
::AxisRotationTransform(unsigned int spaceDimension,
                        unsigned int parametersDimension):
    Superclass(spaceDimension, parametersDimension)
{
    m_Angle = m_TranslationY = m_TranslationZ = itk::NumericTraits<ScalarType>::Zero;
    m_ReferencePlane.SetIdentity();
    m_InverseReferencePlane.SetIdentity();

    this->ComputeMatrix();
}

// Set Parameters
template <class TScalarType>
void
AxisRotationTransform<TScalarType>
::SetParameters(const ParametersType & parameters)
{
    m_Angle = parameters[0];
    m_TranslationY = parameters[1];
    m_TranslationZ = parameters[2];

    this->ComputeMatrix();
    // Modified is always called since we just have a pointer to the
    // parameters and cannot know if the parameters have changed.
    this->Modified();
}


// Get Parameters
template <class TScalarType>
const typename AxisRotationTransform<TScalarType>::ParametersType &
AxisRotationTransform<TScalarType>
::GetParameters() const
{
    this->m_Parameters[0] = m_Angle;
    this->m_Parameters[1] = m_TranslationY;
    this->m_Parameters[2] = m_TranslationZ;

    return this->m_Parameters;
}

// Compose
template <class TScalarType>
void
AxisRotationTransform<TScalarType>
::SetIdentity()
{
    Superclass::SetIdentity();
    m_Angle = m_TranslationY = m_TranslationZ = itk::NumericTraits<ScalarType>::Zero;

    this->ComputeMatrix();
}

// Compute angles from the rotation matrix
template <class TScalarType>
void
AxisRotationTransform<TScalarType>
::ComputeMatrix()
{
    double cosa, sina;

    cosa = std::cos(m_Angle);
    sina = std::sin(m_Angle);

    HomogeneousMatrixType insideMatrix, resultMatrix;
    insideMatrix.SetIdentity();

    insideMatrix[1][1] = cosa;
    insideMatrix[1][2] = -sina;

    insideMatrix[2][1] = sina;
    insideMatrix[2][2] = cosa;

    insideMatrix[1][3] = m_TranslationY;
    insideMatrix[2][3] = m_TranslationZ;

    resultMatrix = m_InverseReferencePlane * insideMatrix * m_ReferencePlane;

    MatrixType linearMatrix;
    OffsetType off;

    for (unsigned int i = 0;i < 3;++i)
    {
        off[i] = resultMatrix(i,3) + this->GetCenter()[i];
        for (unsigned int j = 0;j < 3;++j)
        {
            off[i] -= resultMatrix(i,j) * this->GetCenter()[j];
            linearMatrix(i,j) = resultMatrix(i,j);
        }
    }

    this->SetMatrix(linearMatrix);
    this->SetOffset(off);
}

template<class TScalarType>
void
AxisRotationTransform<TScalarType>::
SetReferencePlaneNormal(InputVectorType &refPlaneNormal)
{
    // The input transforms are matrix transforms realigning each image onto the mid-sagittal plane
    InputVectorType xAxis;
    xAxis.Fill(0.0);
    xAxis[0] = 1;

    MatrixType rotationMatrix = anima::GetRotationMatrixFromVectors(refPlaneNormal,xAxis);
    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j < 3;++j)
            m_ReferencePlane(i,j) = rotationMatrix(i,j);

    m_InverseReferencePlane = m_ReferencePlane.GetInverse();
}

// Print self
template<class TScalarType>
void
AxisRotationTransform<TScalarType>::
PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os,indent);

    os << indent << "Symmetry plane parameters: Angle=" << m_Angle
       << " TranslationY=" << m_TranslationY
       << " TranslationZ=" << m_TranslationZ
       << std::endl;
}

} // end of namespace anima
