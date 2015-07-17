#pragma once

#include <itkMatrix.h>

namespace anima
{

template <class VectorType> itk::Matrix <double,3,3> GetRotationMatrixFromVectors(const VectorType &first_direction, const VectorType &second_direction, const unsigned int dimension);
template <class ScalarType> itk::Matrix <double,3,3> GetRotationMatrixFromVectors(const itk::Point<ScalarType> &first_direction, const itk::Point<ScalarType> &second_direction);
template <class ScalarType, unsigned int NDimension> itk::Matrix <double,3,3> GetRotationMatrixFromVectors(const itk::Vector<ScalarType,NDimension> &first_direction, const itk::Vector<ScalarType,NDimension> &second_direction);
template <class ScalarType, unsigned int NDimension> itk::Matrix <double,3,3> GetRotationMatrixFromVectors(const vnl_vector_fixed<ScalarType,NDimension> &first_direction, const vnl_vector_fixed<ScalarType,NDimension> &second_direction);

} // end of namespace anima

#include "animaMatrixOperations.hxx"
