#pragma once

#include <itkMatrix.h>

namespace anima
{

template <class VectorType> itk::Matrix <double,3,3> GetRotationMatrixFromVectors(const VectorType &first_direction, const VectorType &second_direction, const unsigned int dimension);
template <class ScalarType> itk::Matrix <double,3,3> GetRotationMatrixFromVectors(const itk::Point<ScalarType> &first_direction, const itk::Point<ScalarType> &second_direction);
template <class ScalarType, unsigned int NDimension> itk::Matrix <double,3,3> GetRotationMatrixFromVectors(const itk::Vector<ScalarType,NDimension> &first_direction, const itk::Vector<ScalarType,NDimension> &second_direction);
template <class ScalarType, unsigned int NDimension> itk::Matrix <double,3,3> GetRotationMatrixFromVectors(const vnl_vector_fixed<ScalarType,NDimension> &first_direction, const vnl_vector_fixed<ScalarType,NDimension> &second_direction);

//! Solves lower triangular system maxtrix * result = rhs using Gauss elimination. Can be used with rhs and result pointing to the same address in which case, result will be erased. result has to be the same size as rhs
//! If the matrix is not of full rank, rank can be set to something different than zero, in which case the solver will take only the first rank rows
template <class ScalarType, class VectorType> void LowerTriangularSolver(vnl_matrix <ScalarType> &matrix, VectorType &rhs, VectorType &result, unsigned int rank = 0);

//! Solves upper triangular system maxtrix * result = rhs using Gauss elimination. Can be used with rhs and result pointing to the same address in which case, result will be erased. result has to be the same size as rhs
//! If the matrix is not of full rank, rank can be set to something different than zero, in which case the solver will take only the first rank rows
template <class ScalarType, class VectorType> void UpperTriangularSolver(vnl_matrix <ScalarType> &matrix, VectorType &rhs, VectorType &result, unsigned int rank = 0);

} // end of namespace anima

#include "animaMatrixOperations.hxx"
