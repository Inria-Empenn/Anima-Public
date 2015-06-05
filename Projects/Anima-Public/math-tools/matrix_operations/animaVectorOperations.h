#pragma once

#include <itkVector.h>
#include <itkVariableLengthVector.h>
#include <itkPoint.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_vector_fixed.h>
#include <vector>

namespace anima
{

/******* Main function ComputeEuclideanDistance *******/
// Main
template <class VectorType> double ComputeEuclideanDistance(const VectorType &x1, const VectorType &x2, const unsigned int NDimension);
// For itkVector
template <class ScalarType, unsigned int NDimension> double ComputeEuclideanDistance(const itk::Vector <ScalarType,NDimension> &x1, const itk::Vector <ScalarType,NDimension> &x2);
// For itkVariableLengthVector
template <class ScalarType> double ComputeEuclideanDistance(const itk::VariableLengthVector <ScalarType> &x1, const itk::VariableLengthVector <ScalarType> &x2);
// For itkPoint
template <class ScalarType, unsigned int NDimension> double ComputeEuclideanDistance(const itk::Point <ScalarType,NDimension> &x1, const itk::Point <ScalarType,NDimension> &x2);
// For vnl_vector
template <class ScalarType> double ComputeEuclideanDistance(const vnl_vector <ScalarType> &x1, const vnl_vector <ScalarType> &x2);
// For vnl_vector_fixed
template <class ScalarType, unsigned int NDimension> double ComputeEuclideanDistance(const vnl_vector_fixed <ScalarType,NDimension> &x1, const vnl_vector_fixed <ScalarType,NDimension> &x2);
// For std::vector
template <class ScalarType> double ComputeEuclideanDistance(const std::vector <ScalarType> &x1, const std::vector <ScalarType> &x2);
/******************************************************/

template <class VectorType> double ComputePointToSetDistance(const VectorType &x, const std::vector <VectorType> &s);
template <class VectorType> double ComputeDirectedHausdorffDistance(const std::vector <VectorType> &s1, const std::vector <VectorType> &s2);
template <class VectorType> double ComputeHausdorffDistance(const std::vector <VectorType> &s1, const std::vector <VectorType> &s2);
template <class VectorType> double ComputeModifiedDirectedHausdorffDistance(const std::vector <VectorType> &s1, const std::vector <VectorType> &s2);
template <class VectorType> double ComputeModifiedHausdorffDistance(const std::vector <VectorType> &s1, const std::vector <VectorType> &s2);

/******* Main function ExponentialSum *******/
// Main
template <class VectorType> double ExponentialSum(const VectorType &x, const unsigned int NDimension);
// For itkVector
template <class ScalarType, unsigned int NDimension> double ExponentialSum(const itk::Vector <ScalarType,NDimension> &x);
// For itkVariableLengthVector
template <class ScalarType> double ExponentialSum(const itk::VariableLengthVector <ScalarType> &x);
// For itkPoint
template <class ScalarType, unsigned int NDimension> double ExponentialSum(const itk::Point <ScalarType,NDimension> &x);
// For vnl_vector
template <class ScalarType> double ExponentialSum(const vnl_vector <ScalarType> &x);
// For vnl_vector_fixed
template <class ScalarType, unsigned int NDimension> double ExponentialSum(const vnl_vector_fixed <ScalarType,NDimension> &x);
// For std::vector
template <class ScalarType> double ExponentialSum(const std::vector <ScalarType> &x);
/********************************************/

/******* Main function ComputeScalarProduct *******/
// Main
template <class VectorType> double ComputeScalarProduct(const VectorType &v1, const VectorType &v2, const unsigned int NDimension);
// For itkVector
template <class ScalarType, unsigned int NDimension> double ComputeScalarProduct(const itk::Vector <ScalarType,NDimension> &v1, const itk::Vector <ScalarType,NDimension> &v2);
// For itkVariableLengthVector
template <class ScalarType> double ComputeScalarProduct(const itk::VariableLengthVector <ScalarType> &v1, const itk::VariableLengthVector <ScalarType> &v2);
// For itkPoint
template <class ScalarType, unsigned int NDimension> double ComputeScalarProduct(const itk::Point <ScalarType,NDimension> &v1, const itk::Point <ScalarType,NDimension> &v2);
// For vnl_vector
template <class ScalarType> double ComputeScalarProduct(const vnl_vector <ScalarType> &v1, const vnl_vector <ScalarType> &v2);
// For vnl_vector_fixed
template <class ScalarType, unsigned int NDimension> double ComputeScalarProduct(const vnl_vector_fixed <ScalarType,NDimension> &v1, const vnl_vector_fixed <ScalarType,NDimension> &v2);
// For std::vector
template <class ScalarType> double ComputeScalarProduct(const std::vector <ScalarType> &v1, const std::vector <ScalarType> &v2);
/**************************************************/

/******* Function ComputeCrossProduct *******/
// Main
template <class VectorType> void ComputeCrossProduct(const VectorType &v1, const VectorType &v2, const unsigned int NDimension, VectorType &resVec);
// For itkVector
template <class ScalarType, unsigned int NDimension> void ComputeCrossProduct(const itk::Vector <ScalarType,NDimension> &v1, const itk::Vector <ScalarType,NDimension> &v2, itk::Vector <ScalarType,NDimension> &resVec);
// For itkVariableLengthVector
template <class ScalarType> void ComputeCrossProduct(const itk::VariableLengthVector <ScalarType> &v1, const itk::VariableLengthVector <ScalarType> &v2, itk::VariableLengthVector <ScalarType> &resVec);
// For itkPoint
template <class ScalarType, unsigned int NDimension> void ComputeCrossProduct(const itk::Point <ScalarType,NDimension> &v1, const itk::Point <ScalarType,NDimension> &v2, itk::Point <ScalarType,NDimension> &resVec);
// For vnl_vector
template <class ScalarType> void ComputeCrossProduct(const vnl_vector <ScalarType> &v1, const vnl_vector <ScalarType> &v2, vnl_vector <ScalarType> &resVec);
// For vnl_vector_fixed
template <class ScalarType, unsigned int NDimension> void ComputeCrossProduct(const vnl_vector_fixed <ScalarType,NDimension> &v1, const vnl_vector_fixed <ScalarType,NDimension> &v2, vnl_vector_fixed <ScalarType,NDimension> &resVec);
// For std::vector
template <class ScalarType> void ComputeCrossProduct(const std::vector <ScalarType> &v1, const std::vector <ScalarType> &v2, std::vector <ScalarType> &resVec);
/********************************************/

template <class VectorType> double ComputeNorm(const VectorType &v);

/******* Function Normalize *******/
// Main
template <class VectorType> void Normalize(const VectorType &v, const unsigned int NDimension, VectorType &resVec);
// For itkVector
template <class ScalarType, unsigned int NDimension> void Normalize(const itk::Vector <ScalarType,NDimension> &v, itk::Vector <ScalarType,NDimension> &resVec);
// For itkVariableLengthVector
template <class ScalarType> void Normalize(const itk::VariableLengthVector <ScalarType> &v, itk::VariableLengthVector <ScalarType> &resVec);
// For itkPoint
template <class ScalarType, unsigned int NDimension> void Normalize(const itk::Point <ScalarType,NDimension> &v, itk::Point <ScalarType,NDimension> &resVec);
// For vnl_vector
template <class ScalarType> void Normalize(const vnl_vector <ScalarType> &v, vnl_vector <ScalarType> &resVec);
// For vnl_vector_fixed
template <class ScalarType, unsigned int NDimension> void Normalize(const vnl_vector_fixed <ScalarType, NDimension> &v, vnl_vector_fixed <ScalarType, NDimension> &resVec);
// For std::vector
template <class ScalarType> void Normalize(const std::vector <ScalarType> &v, std::vector <ScalarType> &resVec);
/**********************************/

template <class T1, class T2> double ComputeAngle(const std::vector <T1> &v1, const std::vector <T2> &v2);
template <class T1, class T2> double ComputeOrientationAngle(const std::vector <T1> &v1, const std::vector <T2> &v2);
template <class ScalarType, unsigned int NDimension>
double
ComputeOrientationAngle(const itk::Vector <ScalarType,NDimension> &v1, const itk::Vector <ScalarType,NDimension> &v2);

template <class VectorType> void Revert(const VectorType &v, const unsigned int vSize, VectorType &resVec);
template <class ScalarType> void Revert(const std::vector <ScalarType> &v, std::vector<ScalarType> &resVec);
template <class ScalarType> void Revert(const itk::Point <ScalarType> &v, itk::Point <ScalarType> &resVec);
template <class ScalarType, unsigned int NDimension> void Revert(const vnl_vector_fixed <ScalarType, NDimension> &v,
                                                                 vnl_vector_fixed <ScalarType, NDimension> &resVec);

template <class VectorType> void TransformCartesianToSphericalCoordinates(const VectorType &v, VectorType &resVec);
template <class VectorType> void TransformSphericalToCartesianCoordinates(const VectorType &v, VectorType &resVec);
template <class T1, class T2, class T3> void ProjectOnOrthogonalPlane(const std::vector <T1> &v, const std::vector <T2> &normalVec, std::vector <T3> &resVec);

template <class T1, class T2, class T3, class T4> void RotateAroundAxis(const std::vector <T1> &v, const T2 &phi, const std::vector <T3> &normalVec, std::vector <T4> &resVec);
template <class Vector3DType, class ScalarType> void RotateAroundAxis(const Vector3DType &v, const ScalarType &phi, const Vector3DType &normalVec, Vector3DType &resVec);

} // end of namespace anima

#include "animaVectorOperations.hxx"
