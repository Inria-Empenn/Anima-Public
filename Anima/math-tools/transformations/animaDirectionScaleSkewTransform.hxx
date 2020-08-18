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
    // Set angles with parameters
    m_LogScale = parameters[0];
    m_LogFirstSkew = parameters[1];
    m_LogSecondSkew = parameters[2];
    m_UniqueTranslation = parameters[3];

    this->ComputeMatrix();
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

    m_LogTransform.fill(itk::NumericTraits<TScalarType>::Zero);
    m_LogVector.Fill(itk::NumericTraits<TScalarType>::Zero);
}

template <class TScalarType>
void
DirectionScaleSkewTransform<TScalarType>
::SetGeometry(HomogeneousMatrixType &matrix, bool update)
{
    m_Geometry = matrix;
    m_GeometryInv = m_Geometry.GetInverse();

    if (update)
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
    // Start by updating log-representations
    this->GetInternalLogarithm();
    this->GetInternalExponential();

    MatrixType linearMatrix;
    OffsetType off;

    // Apply geometry to log and exp transforms
    this->InternalApplyGeometry(linearMatrix,off);

    this->SetMatrix(linearMatrix);
    this->SetOffset(off);

    m_LogVector.Fill(itk::NumericTraits<TScalarType>::Zero);
    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j <= 3;++j)
            m_LogVector[i * 4 + j] = m_LogTransform(i,j);
}

template <class TScalarType>
void
DirectionScaleSkewTransform<TScalarType>
::InternalApplyGeometry(MatrixType &linearMatrix, OffsetType &offset)
{
    // First update log transform
    for (unsigned int i = 0;i < 3;++i)
        offset[i] = m_Geometry(i,m_UniqueDirection) * m_UniqueTranslation;

    for (unsigned int j = 0;j < 3;++j)
    {
        double aRinv = 0.0;
        for (unsigned int k = 0;k < 3;++k)
            aRinv += m_LogTransform(m_UniqueDirection,k) * m_GeometryInv(k,j);

        for (unsigned int i = 0;i < 3;++i)
            linearMatrix(i,j) = aRinv * m_Geometry(i,m_UniqueDirection);
    }

    for (unsigned int i = 0;i < 3;++i)
    {
        for (unsigned int j = 0;j < 3;++j)
            offset[i] -= linearMatrix(i,j) * m_Geometry(j,3);

        m_LogTransform(i,3) = offset[i];
        for (unsigned int j = 0;j < 3;++j)
            m_LogTransform(i,j) = linearMatrix(i,j);
    }

    // Then update exp transform and keep linearMatrix and offset there
    for (unsigned int i = 0;i < 3;++i)
        offset[i] = m_Geometry(i,m_UniqueDirection) * m_ExpTransform(m_UniqueDirection,3) + m_Geometry(i,3);

    for (unsigned int j = 0;j < 3;++j)
    {
        double aRinv = 0.0;
        for (unsigned int k = 0;k < 3;++k)
            aRinv += m_ExpTransform(m_UniqueDirection,k) * m_GeometryInv(k,j);

        for (unsigned int i = 0;i < 3;++i)
        {
            linearMatrix(i,j) = aRinv * m_Geometry(i,m_UniqueDirection);
            for (unsigned int k = 0;k < 3;++k)
            {
                if (k == m_UniqueDirection)
                    continue;
                linearMatrix(i,j) += m_Geometry(i,k) * m_GeometryInv(k,j);
            }
        }
    }

    for (unsigned int i = 0;i < 3;++i)
    {
        for (unsigned int j = 0;j < 3;++j)
            offset[i] -= linearMatrix(i,j) * m_Geometry(j,3);

        m_ExpTransform(i,3) = offset[i];
        for (unsigned int j = 0;j < 3;++j)
            m_ExpTransform(i,j) = linearMatrix(i,j);
    }
}

template <class TScalarType>
void
DirectionScaleSkewTransform<TScalarType>
::GetInternalLogarithm()
{
    m_LogTransform.fill(0.0);

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
}

template <class TScalarType>
void
DirectionScaleSkewTransform<TScalarType>
::GetInternalExponential()
{
    m_ExpTransform.set_identity();

    unsigned int skewLocation = 0;
    if (m_UniqueDirection == skewLocation)
        skewLocation++;

    double expFactor = 1.0;
    if (m_LogScale != 0.0)
    {
        expFactor = std::exp(m_LogScale);
        m_ExpTransform(m_UniqueDirection,skewLocation) = m_LogFirstSkew * (expFactor - 1.0) / m_LogScale;
    }
    else
        m_ExpTransform(m_UniqueDirection,skewLocation) = m_LogFirstSkew;

    skewLocation++;
    if (m_UniqueDirection == skewLocation)
        skewLocation++;

    if (m_LogScale != 0.0)
        m_ExpTransform(m_UniqueDirection,skewLocation) = m_LogSecondSkew * (expFactor - 1.0) / m_LogScale;
    else
        m_ExpTransform(m_UniqueDirection,skewLocation) = m_LogSecondSkew;

    m_ExpTransform(m_UniqueDirection,m_UniqueDirection) = expFactor;

    if (m_LogScale != 0.0)
        m_ExpTransform(m_UniqueDirection,3) = m_UniqueTranslation * (expFactor - 1.0) / m_LogScale;
    else
        m_ExpTransform(m_UniqueDirection,3) = m_UniqueTranslation;
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

//// Direction scale transform from here

// Set Parameters
template <class TScalarType>
void
DirectionScaleTransform<TScalarType>
::SetParameters(const ParametersType & parameters)
{
    // Set angles with parameters
    this->SetLogScale(parameters[0],false);
    this->SetLogFirstSkew(0,false);
    this->SetLogSecondSkew(0,false);
    this->SetUniqueTranslation(parameters[1],false);

    this->ComputeMatrix();
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

template <class TScalarType>
void
DirectionScaleTransform<TScalarType>
::InternalApplyGeometry(MatrixType &linearMatrix, OffsetType &offset)
{
    unsigned int uniqueDirection = this->GetUniqueDirection();
    // First update log transform
    for (unsigned int i = 0;i < 3;++i)
        offset[i] = this->m_Geometry(i,uniqueDirection) * this->GetUniqueTranslation();

    for (unsigned int j = 0;j < 3;++j)
    {
        double aRinv = this->m_LogTransform(uniqueDirection,uniqueDirection) * this->m_GeometryInv(uniqueDirection,j);

        for (unsigned int i = 0;i < 3;++i)
            linearMatrix(i,j) = aRinv * this->m_Geometry(i,uniqueDirection);
    }

    for (unsigned int i = 0;i < 3;++i)
    {
        for (unsigned int j = 0;j < 3;++j)
            offset[i] -= linearMatrix(i,j) * this->m_Geometry(j,3);

        this->m_LogTransform(i,3) = offset[i];
        for (unsigned int j = 0;j < 3;++j)
            this->m_LogTransform(i,j) = linearMatrix(i,j);
    }

    // Then update exp transform and keep linearMatrix and offset there
    for (unsigned int i = 0;i < 3;++i)
        offset[i] = this->m_Geometry(i,uniqueDirection) * this->m_ExpTransform(uniqueDirection,3) + this->m_Geometry(i,3);

    for (unsigned int j = 0;j < 3;++j)
    {
        double aRinv = this->m_ExpTransform(uniqueDirection,uniqueDirection) * this->m_GeometryInv(uniqueDirection,j);

        for (unsigned int i = 0;i < 3;++i)
        {
            linearMatrix(i,j) = aRinv * this->m_Geometry(i,uniqueDirection);
            for (unsigned int k = 0;k < 3;++k)
            {
                if (k == uniqueDirection)
                    continue;
                linearMatrix(i,j) += this->m_Geometry(i,k) * this->m_GeometryInv(k,j);
            }
        }
    }

    for (unsigned int i = 0;i < 3;++i)
    {
        for (unsigned int j = 0;j < 3;++j)
            offset[i] -= linearMatrix(i,j) * this->m_Geometry(j,3);

        this->m_ExpTransform(i,3) = offset[i];
        for (unsigned int j = 0;j < 3;++j)
            this->m_ExpTransform(i,j) = linearMatrix(i,j);
    }
}

template <class TScalarType>
void
DirectionScaleTransform<TScalarType>
::GetInternalLogarithm()
{
    this->m_LogTransform.fill(0.0);

    unsigned int uniqueDirection = this->GetUniqueDirection();
    this->m_LogTransform(uniqueDirection,uniqueDirection) = this->GetLogScale();
    this->m_LogTransform(uniqueDirection,3) = this->GetUniqueTranslation();
}

template <class TScalarType>
void
DirectionScaleTransform<TScalarType>
::GetInternalExponential()
{
    this->m_ExpTransform.set_identity();

    double expFactor = 1.0;
    double logScale = this->GetLogScale();
    if (logScale != 0.0)
        expFactor = std::exp(logScale);

    unsigned int uniqueDirection = this->GetUniqueDirection();
    this->m_ExpTransform(uniqueDirection,uniqueDirection) = expFactor;

    if (logScale != 0.0)
        this->m_ExpTransform(uniqueDirection,3) = this->GetUniqueTranslation() * (expFactor - 1.0) / logScale;
    else
        this->m_ExpTransform(uniqueDirection,3) = this->GetUniqueTranslation();
}

//// Direction transform from here

// Set Parameters
template <class TScalarType>
void
DirectionTransform<TScalarType>
::SetParameters(const ParametersType & parameters)
{
    // Set angles with parameters
    this->SetLogScale(0,false);
    this->SetLogFirstSkew(0,false);
    this->SetLogSecondSkew(0,false);
    this->SetUniqueTranslation(parameters[0],false);

    this->ComputeMatrix();
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

template <class TScalarType>
void
DirectionTransform<TScalarType>
::InternalApplyGeometry(MatrixType &linearMatrix, OffsetType &offset)
{
    unsigned int uniqueDirection = this->GetUniqueDirection();
    // First update log transform
    for (unsigned int i = 0;i < 3;++i)
        offset[i] = this->m_Geometry(i,uniqueDirection) * this->GetUniqueTranslation();

    for (unsigned int i = 0;i < 3;++i)
        this->m_LogTransform(i,3) = offset[i];

    // Then update exp transform and keep linearMatrix and offset there
    for (unsigned int i = 0;i < 3;++i)
        offset[i] = this->m_Geometry(i,uniqueDirection) * this->m_ExpTransform(uniqueDirection,3);

    linearMatrix.SetIdentity();

    for (unsigned int i = 0;i < 3;++i)
    {
        this->m_ExpTransform(i,3) = offset[i];
        for (unsigned int j = 0;j < 3;++j)
            this->m_ExpTransform(i,j) = linearMatrix(i,j);
    }
}

template <class TScalarType>
void
DirectionTransform<TScalarType>
::GetInternalLogarithm()
{
    this->m_LogTransform.fill(0.0);
    this->m_LogTransform(this->GetUniqueDirection(),3) = this->GetUniqueTranslation();
}

template <class TScalarType>
void
DirectionTransform<TScalarType>
::GetInternalExponential()
{
    this->m_ExpTransform.set_identity();
    this->m_ExpTransform(this->GetUniqueDirection(),3) = this->GetUniqueTranslation();
}

} // end namespace anima
