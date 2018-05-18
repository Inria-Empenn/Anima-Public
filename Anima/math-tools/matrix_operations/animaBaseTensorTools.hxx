#pragma once

#include "animaBaseTensorTools.h"
#include <itkSymmetricEigenAnalysis.h>

#include <animaVectorOperations.h>
#include <animaMatrixOperations.h>

namespace anima
{
template <class T>
void
GetTensorLogarithm(const vnl_matrix <T> &tensor, vnl_matrix <T> &log_tensor)
{
    typedef itk::SymmetricEigenAnalysis < vnl_matrix <T>, vnl_diag_matrix<T>, vnl_matrix <T> > EigenAnalysisType;
    unsigned int tensDim = tensor.rows();

    EigenAnalysisType eigen(tensDim);
    vnl_matrix <T> eigVecs(tensDim,tensDim);
    vnl_diag_matrix <T> eigVals(tensDim);

    eigen.ComputeEigenValuesAndVectors(tensor,eigVals,eigVecs);

    for (unsigned int i = 0;i < tensDim;++i)
    {
        if (eigVals[i] <= 1.0e-16)
            eigVals[i] = 1.0e-16;

        eigVals[i] = std::log(eigVals[i]);
    }

    RecomposeTensor(eigVals,eigVecs,log_tensor);
}

template <class T>
void
GetTensorPower(const vnl_matrix <T> &tensor, vnl_matrix <T> &outputTensor, double powerValue)
{
    typedef itk::SymmetricEigenAnalysis < vnl_matrix <T>, vnl_diag_matrix<T>, vnl_matrix <T> > EigenAnalysisType;
    unsigned int tensDim = tensor.rows();

    EigenAnalysisType eigen(tensDim);
    vnl_matrix <T> eigVecs(tensDim,tensDim);
    vnl_diag_matrix <T> eigVals(tensDim);

    eigen.ComputeEigenValuesAndVectors(tensor,eigVals,eigVecs);

    for (unsigned int i = 0;i < tensDim;++i)
    {
        if (eigVals[i] <= 1.0e-16)
            eigVals[i] = 1.0e-16;

        eigVals[i] = std::pow(eigVals[i],powerValue);
    }

    RecomposeTensor(eigVals,eigVecs,outputTensor);
}

template <class T>
void
GetTensorExponential(const vnl_matrix <T> &log_tensor, vnl_matrix <T> &tensor)
{
    typedef itk::SymmetricEigenAnalysis < vnl_matrix <T>, vnl_diag_matrix<T>, vnl_matrix <T> > EigenAnalysisType;
    unsigned int tensDim = log_tensor.rows();

    EigenAnalysisType eigen(tensDim);
    vnl_matrix <T> eigVecs(tensDim,tensDim);
    vnl_diag_matrix <T> eigVals(tensDim);

    eigen.ComputeEigenValuesAndVectors(log_tensor,eigVals,eigVecs);

    for (unsigned int i = 0;i < tensDim;++i)
        eigVals[i] = std::exp(eigVals[i]);

    RecomposeTensor(eigVals,eigVecs,tensor);
}

template <class T1, class T2>
void
GetVectorRepresentation(const vnl_matrix <T1> &tensor, itk::VariableLengthVector <T2> &vector,
                        unsigned int vecDim, bool scale)
{
    unsigned int dim = tensor.rows();
    if (vecDim == 0)
        vecDim = dim * (dim + 1) / 2;

    double sqrt2 = sqrt(2.0);
    if (vector.GetSize() != vecDim)
        vector.SetSize(vecDim);

    unsigned int pos = 0;
    for (unsigned int i = 0;i < dim;++i)
        for (unsigned int j = 0;j <= i;++j)
        {
            vector[pos] = tensor(i,j);
            if ((i != j)&&scale)
                vector[pos] *= sqrt2;
            ++pos;
        }
}

template <class T1, class T2>
void
GetTensorFromVectorRepresentation(const itk::VariableLengthVector <T1> &vector,
                                  vnl_matrix <T2> &tensor, unsigned int tensDim, bool scale)
{
    unsigned int vecDim = vector.GetSize();

    if (tensDim == 0)
        tensDim = std::floor((std::sqrt((float)(8 * vecDim + 1)) - 1) / 2.0);

    double sqrt2 = std::sqrt(2.0);
    tensor.set_size(tensDim,tensDim);

    unsigned int pos = 0;
    for (unsigned int i = 0;i < tensDim;++i)
        for (unsigned int j = 0;j <= i;++j)
        {
            tensor(i,j) = vector[pos];
            if (i != j)
            {
                if (scale)
                    tensor(i,j) /= sqrt2;

                tensor(j,i) = tensor(i,j);
            }

            ++pos;
        }
}

template <class T> void ProjectOnTensorSpace(const vnl_matrix <T> &matrix, vnl_matrix <T> &tensor)
{
    typedef itk::SymmetricEigenAnalysis < vnl_matrix <T>, vnl_diag_matrix<T>, vnl_matrix <T> > EigenAnalysisType;
    unsigned int tensDim = matrix.rows();

    EigenAnalysisType eigen(tensDim);
    vnl_matrix <T> eigVecs(tensDim,tensDim);
    vnl_diag_matrix <T> eigVals(tensDim);

    eigen.ComputeEigenValuesAndVectors(matrix,eigVals,eigVecs);

    for (unsigned int i = 0;i < tensDim;++i)
    {
        if (eigVals[i] <= 1.0e-16)
            eigVals[i] = 1.0e-16;
    }

    RecomposeTensor(eigVals,eigVecs,tensor);
}

template <typename RealType>
void
ExtractRotationFromJacobianMatrix(vnl_matrix <RealType> &jacobianMatrix, vnl_matrix <RealType> &rotationMatrix,
                                  vnl_matrix <RealType> &tmpMat)
{
    unsigned int tensorDimension = 3;
    tmpMat.set_size(tensorDimension, tensorDimension);

    for (unsigned int l = 0;l < tensorDimension;++l)
        for (unsigned int m = l;m < tensorDimension;++m)
        {
            tmpMat(l,m) = 0;
            for (unsigned int n = 0;n < tensorDimension;++n)
                tmpMat(l,m) += jacobianMatrix(l,n)*jacobianMatrix(m,n);

            if (m != l)
                tmpMat(m,l) = tmpMat(l,m);
        }

    itk::SymmetricEigenAnalysis < vnl_matrix <RealType>, vnl_diag_matrix<RealType>, vnl_matrix <RealType> > eigenComputer(tensorDimension);
    vnl_matrix <RealType> eVec(tensorDimension,tensorDimension);
    vnl_diag_matrix <RealType> eVals(tensorDimension);

    eigenComputer.ComputeEigenValuesAndVectors(tmpMat, eVals, eVec);

    for (unsigned int i = 0;i < tensorDimension;++i)
        eVals[i] = pow(eVals[i], -0.5);

    RecomposeTensor(eVals,eVec,rotationMatrix);

    rotationMatrix *= jacobianMatrix;
}

template <typename RealType, typename MatrixType>
void
ExtractPPDRotationFromJacobianMatrix(vnl_matrix <RealType> &jacobianMatrix, vnl_matrix <RealType> &rotationMatrix,
                                     MatrixType &eigenVectors)
{
    const unsigned int tensorDimension = 3;

    typedef vnl_vector_fixed <RealType,tensorDimension> VectorType;
    VectorType transformedPrincipal;
    VectorType transformedSecondary;
    for (unsigned int i = 0;i < tensorDimension;++i)
    {
        transformedPrincipal[i] = 0;
        transformedSecondary[i] = 0;
        for (unsigned int j = 0;j < tensorDimension;++j)
        {
            transformedPrincipal[i] += jacobianMatrix(i,j) * eigenVectors(2,j);
            transformedSecondary[i] += jacobianMatrix(i,j) * eigenVectors(1,j);
        }
    }

    anima::Normalize(transformedPrincipal,transformedPrincipal);
    anima::Normalize(transformedSecondary,transformedSecondary);

    VectorType principalEigenVector;
    for (unsigned int i = 0;i < tensorDimension;++i)
        principalEigenVector[i] = eigenVectors(2,i);
    MatrixType firstRotation = anima::GetRotationMatrixFromVectors(principalEigenVector,transformedPrincipal);

    VectorType rotatedSecondary;
    double dotProductTransformed = 0;
    for (unsigned int i = 0;i < tensorDimension;++i)
    {
        dotProductTransformed += transformedPrincipal[i] * transformedSecondary[i];
        rotatedSecondary[i] = 0;
        for (unsigned int j = 0;j < tensorDimension;++j)
            rotatedSecondary[i] += firstRotation(i,j) * eigenVectors(1,j);
    }

    VectorType projectedSecondary;
    for (unsigned int i = 0;i < tensorDimension;++i)
        projectedSecondary[i] = transformedSecondary[i] - dotProductTransformed * transformedPrincipal[i];

    anima::Normalize(projectedSecondary,projectedSecondary);

    MatrixType secondRotation = anima::GetRotationMatrixFromVectors(rotatedSecondary,projectedSecondary);
    rotationMatrix.set_size(tensorDimension,tensorDimension);
    rotationMatrix.fill(0.0);
    for (unsigned int i = 0;i < tensorDimension;++i)
    {
        for (unsigned int j = 0;j < tensorDimension;++j)
        {
            for (unsigned int k = 0;k < tensorDimension;++k)
                rotationMatrix(i,j) += secondRotation(i,k) * firstRotation(k,j);
        }
    }
}

template <class T1, class T2>
void
RecomposeTensor(vnl_diag_matrix <T1> &eigs, vnl_matrix <T1> &eigVecs, vnl_matrix <T2> &resMatrix)
{
    unsigned int tensorDimension = eigs.size();
    resMatrix.set_size(tensorDimension,tensorDimension);

    for (unsigned int i = 0;i < tensorDimension;++i)
    {
        for (unsigned int j = i;j < tensorDimension;++j)
        {
            resMatrix(i,j) = 0;
            for (unsigned int k = 0;k < tensorDimension;++k)
                resMatrix(i,j) += eigs[k] * eigVecs(k,i) * eigVecs(k,j);

            if (j != i)
                resMatrix(j,i) = resMatrix(i,j);
        }
    }
}

template <class T1, class T2, class T3>
void RotateSymmetricMatrix(T1 &tensor, T2 &rotationMatrix, T3 &rotated_tensor, unsigned int tensorDim)
{
    for (unsigned int i = 0;i < tensorDim;++i)
        for (unsigned int j = i;j < tensorDim;++j)
        {
            rotated_tensor(i,j) = 0;
            for (unsigned int a = 0;a < tensorDim;++a)
                for (unsigned int b = a;b < tensorDim;++b)
                {
                    double factor = rotationMatrix(b,i) * rotationMatrix(a,j);
                    if (b != a)
                        factor += rotationMatrix(a,i) * rotationMatrix(b,j);

                    rotated_tensor(i,j) += tensor(a,b) * factor;
                }

            if (i != j)
                rotated_tensor(j,i) = rotated_tensor(i,j);
        }
}

template <class T1, class T2>
void RotateSymmetricMatrix(vnl_matrix <T1> &tensor, vnl_matrix <T2> &rotationMatrix, vnl_matrix <T2> &rotated_tensor)
{
    unsigned int tensorDim = tensor.rows();
    rotated_tensor.set_size(tensorDim,tensorDim);

    anima::RotateSymmetricMatrix(tensor,rotationMatrix,rotated_tensor,tensorDim);
}

template <class T1, class T2, unsigned int NDim>
void RotateSymmetricMatrix(itk::Matrix <T1,NDim,NDim> &tensor, itk::Matrix <T2,NDim,NDim> &rotationMatrix,
                           itk::Matrix <T2,NDim,NDim> &rotated_tensor)
{
    anima::RotateSymmetricMatrix(tensor,rotationMatrix,rotated_tensor,NDim);
}

template <class T1>
double
ovlScore(vnl_diag_matrix <T1> &eigsX, vnl_matrix <T1> &eigVecsX,
         vnl_diag_matrix <T1> &eigsY, vnl_matrix <T1> &eigVecsY)
{
    unsigned int sizeData = eigsX.size();

    long double resVal = 0;
    long double denomVal = 0;

    for (unsigned int i = 0;i < sizeData;++i)
    {
        long double dotProd = 0;
        for (unsigned int j = 0;j < sizeData;++j)
            dotProd += eigVecsX(i,j)*eigVecsY(i,j);

        resVal += dotProd * dotProd * eigsX[i] * eigsY[i];
        denomVal += eigsX[i] * eigsY[i];
    }

    return resVal/denomVal;
}


} // end of namespace anima
