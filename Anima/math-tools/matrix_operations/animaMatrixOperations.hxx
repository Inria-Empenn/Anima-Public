#pragma once
#include "animaMatrixOperations.h"

#include <animaVectorOperations.h>
#include <animaLinearTransformEstimationTools.h>
#include <itkSymmetricEigenAnalysis.h>

namespace anima
{

template <class VectorType>
itk::Matrix <double,3,3>
GetRotationMatrixFromVectors(const VectorType &first_direction, const VectorType &second_direction, const unsigned int dimension)
{
    vnl_matrix<double> tmpMatrix(dimension+1,dimension+1,0);
    anima::pairingToQuaternion(first_direction,second_direction,tmpMatrix);
    vnl_matrix<double> AMatrix = tmpMatrix.transpose()*tmpMatrix;

    // Needs two points, providing one on normal vector (cross product)
    VectorType tmpVec;
    tmpVec[0] = first_direction[1] * second_direction[2] - first_direction[2] * second_direction[1];
    tmpVec[1] = first_direction[2] * second_direction[0] - first_direction[0] * second_direction[2];
    tmpVec[2] = first_direction[0] * second_direction[1] - first_direction[1] * second_direction[0];

    double normTmpVec = 0;
    for (unsigned int i = 0;i < dimension;++i)
        normTmpVec += tmpVec[i] * tmpVec[i];

    normTmpVec = std::sqrt(normTmpVec);
    if (normTmpVec < 1.0e-8)
    {
        tmpVec[0] = 0;
        tmpVec[1] = 0;
        tmpVec[2] = 1;
    }

    anima::pairingToQuaternion(tmpVec,tmpVec,tmpMatrix);
    AMatrix += tmpMatrix.transpose()*tmpMatrix;

    itk::SymmetricEigenAnalysis < vnl_matrix <double>, vnl_diag_matrix<double>, vnl_matrix <double> > eigenSystem(dimension+1);
    vnl_matrix <double> eVec(dimension+1,dimension+1);
    vnl_diag_matrix <double> eVals(dimension+1);

    eigenSystem.SetOrderEigenValues(true);
    eigenSystem.ComputeEigenValuesAndVectors(AMatrix, eVals, eVec);

    vnl_matrix <double> rotationMatrix = anima::computeRotationFromQuaternion<double,double>(eVec.get_row(0));

    return rotationMatrix;
}

template <class ScalarType>
itk::Matrix <double,3,3>
GetRotationMatrixFromVectors(const itk::Point<ScalarType> &first_direction, const itk::Point<ScalarType> &second_direction)
{
    unsigned int dimension = first_direction.GetPointDimension();
    return GetRotationMatrixFromVectors(first_direction, second_direction, dimension);
}

template <class ScalarType, unsigned int NDimension>
itk::Matrix <double,3,3>
GetRotationMatrixFromVectors(const itk::Vector<ScalarType,NDimension> &first_direction, const itk::Vector<ScalarType,NDimension> &second_direction)
{
    return GetRotationMatrixFromVectors(first_direction, second_direction, NDimension);
}

template <class ScalarType, unsigned int NDimension>
itk::Matrix <double,3,3>
GetRotationMatrixFromVectors(const vnl_vector_fixed<ScalarType,NDimension> &first_direction, const vnl_vector_fixed<ScalarType,NDimension> &second_direction)
{
    return GetRotationMatrixFromVectors(first_direction, second_direction, NDimension);
}

} // end of namespace anima
