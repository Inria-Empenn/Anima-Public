#pragma once

#include <itkMatrix.h>

namespace anima
{

typedef itk::Matrix <double,3,3> MatrixType;

template <class VectorType> MatrixType GetRotationMatrixFromVectors(const VectorType &first_direction, const VectorType &second_direction, const unsigned int dimension);
template <class ScalarType> MatrixType GetRotationMatrixFromVectors(const std::vector <ScalarType> &first_direction, const std::vector <ScalarType> &second_direction);
template <class ScalarType> MatrixType GetRotationMatrixFromVectors(const itk::Point<ScalarType> &first_direction, const itk::Point<ScalarType> &second_direction);
template <class ScalarType, unsigned int NDimension> MatrixType GetRotationMatrixFromVectors(const itk::Vector<ScalarType,NDimension> &first_direction, const itk::Vector<ScalarType,NDimension> &second_direction);

} // end of namespace anima

#include "animaMatrixOperations.hxx"
