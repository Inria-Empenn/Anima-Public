#pragma once

#include "animaBaseTensorTools.h"
#include <itkSymmetricEigenAnalysis.h>

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
            eigVals[i] = log(1.0e-16);
        else
            eigVals[i] = log(eigVals[i]);
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
        eigVals[i] = std::pow(eigVals[i],powerValue);

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
        eigVals[i] = exp(eigVals[i]);

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
        tensDim = floor((sqrt((float)(8 * vecDim + 1)) - 1) / 2.0);

    double sqrt2 = sqrt(2.0);
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

template <typename PixelRawType, unsigned int ImageDimension>
void
ExtractRotationFromMatrixTransform(itk::MatrixOffsetTransformBase <PixelRawType,ImageDimension,ImageDimension> *trsf,
                                   vnl_matrix <double> &rotationMatrix, vnl_matrix <double> &tmpMat)
{
    unsigned int tensorDimension = 3;

    vnl_matrix <double> jacMatrix(tensorDimension, tensorDimension);
    tmpMat.set_size(tensorDimension, tensorDimension);

    for (unsigned int i = 0;i < tensorDimension;++i)
        for (unsigned int j = 0;j < tensorDimension;++j)
            jacMatrix(i,j) = trsf->GetMatrix()(i,j);

    for (unsigned int l = 0;l < tensorDimension;++l)
        for (unsigned int m = l;m < tensorDimension;++m)
        {
            tmpMat(l,m) = 0;
            for (unsigned int n = 0;n < tensorDimension;++n)
                tmpMat(l,m) += jacMatrix(l,n)*jacMatrix(m,n);

            if (m != l)
                tmpMat(m,l) = tmpMat(l,m);
        }

    itk::SymmetricEigenAnalysis < vnl_matrix <double>, vnl_diag_matrix<double>, vnl_matrix <double> > eigenComputer(tensorDimension);
    vnl_matrix <double> eVec(tensorDimension,tensorDimension);
    vnl_diag_matrix <double> eVals(tensorDimension);

    eigenComputer.ComputeEigenValuesAndVectors(tmpMat, eVals, eVec);

    for (unsigned int i = 0;i < tensorDimension;++i)
        eVals[i] = pow(eVals[i], -0.5);

    RecomposeTensor(eVals,eVec,rotationMatrix);

    rotationMatrix *= jacMatrix;
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

template <class T1, class T2>
void
RotateSymmetricMatrix(vnl_matrix <T1> &tensor, vnl_matrix <T2> &rotationMatrix, vnl_matrix <T2> &rotated_tensor)
{
    //tensor = rotationMatrixTranspose * tensor;
    //tensor *= rotationMatrix;
    unsigned int tensorDim = tensor.rows();
    rotated_tensor.set_size(tensorDim,tensorDim);

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
