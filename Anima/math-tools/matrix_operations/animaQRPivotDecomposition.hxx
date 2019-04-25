#pragma once
#include "animaQRPivotDecomposition.h"
#include <animaVectorOperations.h>

namespace anima
{

template <typename ScalarType> void QRPivotDecomposition(vnl_matrix <ScalarType> &aMatrix, std::vector <unsigned int> &pivotVector,
                                                         std::vector <ScalarType> &houseBetaValues, unsigned int &rank)
{
    unsigned int m = aMatrix.rows();
    unsigned int n = aMatrix.cols();

    if (m < n)
        return;

    pivotVector.resize(n);
    houseBetaValues.resize(n);

    for (unsigned int i = 0;i < n;++i)
    {
        pivotVector[i] = i;
        houseBetaValues[i] = 0.0;
    }

    std::vector <double> cVector(n,0.0);
    for (unsigned int i = 0;i < n;++i)
    {
        for (unsigned int j = 0;j < m;++j)
            cVector[i] += aMatrix(j,i) * aMatrix(j,i);
    }

    rank = 0;
    double tau = cVector[0];
    unsigned int k = 0;
    for (unsigned int i = 1;i < n;++i)
    {
        if (tau < cVector[i])
        {
            k = i;
            tau = cVector[i];
        }
    }

    std::vector <double> vVector;
    std::vector <double> housedVector;
    vnl_matrix <double> workMatrixHouse(m,m);
    vnl_matrix <double> workMatrix(m,n);

    while (tau > 1.0e-16)
    {
        ++rank;
        unsigned int r = rank - 1;
        unsigned int tmpIndex = pivotVector[k];
        pivotVector[k] = pivotVector[r];
        pivotVector[r] = tmpIndex;
        double tmp = cVector[k];
        cVector[k] = cVector[r];
        cVector[r] = tmp;

        for (unsigned int i = 0;i < m;++i)
        {
            tmp = aMatrix(i,k);
            aMatrix(i,k) = aMatrix(i,r);
            aMatrix(i,r) = tmp;
        }

        housedVector.resize(m - r);
        for (unsigned int i = r;i < m;++i)
            housedVector[i - r] = aMatrix(i,r);

        ComputeHouseholderVector(housedVector,vVector,houseBetaValues[r]);

        for (unsigned int i = r;i < m;++i)
        {
            for (unsigned int j = r;j < n;++j)
                workMatrix(i,j) = aMatrix(i,j);
        }

        for (unsigned int i = 0;i < m - r;++i)
        {
            for (unsigned int j = i;j < m - r;++j)
            {
                workMatrixHouse(i,j) = - houseBetaValues[r] * vVector[i] * vVector[j];
                if (i != j)
                    workMatrixHouse(j,i) = workMatrixHouse(i,j);
                else
                    workMatrixHouse(i,j) += 1.0;
            }
        }

        for (unsigned int i = r;i < m;++i)
        {
            unsigned int iIndex = i - r;
            for (unsigned int j = r;j < n;++j)
            {
                aMatrix(i,j) = 0.0;
                for (unsigned int l = r;l < m;++l)
                {
                    unsigned int lIndex = l - r;
                    aMatrix(i,j) += workMatrixHouse(iIndex,lIndex) * workMatrix(l,j);
                }
            }
        }

        for (unsigned int i = rank;i < m;++i)
            aMatrix(i,r) = vVector[i - r];

        for (unsigned int i = rank;i < n;++i)
            cVector[i] -= aMatrix(r,i) * aMatrix(r,i);

        tau = 0.0;

        if (rank < n)
        {
            k = rank;
            tau = cVector[rank];

            for (unsigned int i = rank + 1;i < n;++i)
            {
                if (tau < cVector[i] - 1.0e-16)
                {
                    k = i;
                    tau = cVector[i];
                }
            }
        }
    }
}

template <typename ScalarType> void GetQtBFromQRDecomposition(vnl_matrix <ScalarType> &qrMatrix, vnl_vector<ScalarType> &bVector,
                                                              std::vector <ScalarType> &houseBetaValues, unsigned int rank)
{
    unsigned int m = qrMatrix.rows();
    unsigned int n = qrMatrix.cols();

    if (m < n)
        return;

    for (unsigned int j = 0;j < rank;++j)
    {
        // Compute v^T B
        double vtB = bVector[j];
        for (unsigned int i = j + 1;i < m;++i)
            vtB += qrMatrix(i,j) * bVector[i];

        // bVector -= beta * v * vtB
        bVector[j] -= houseBetaValues[j] * vtB;

        for (unsigned int i = j + 1;i < m;++i)
            bVector[i] -= houseBetaValues[j] * qrMatrix(i,j) * vtB;
    }
}

template <typename ScalarType> void GetQMatrixQRDecomposition(vnl_matrix <ScalarType> &qrMatrix, std::vector <ScalarType> &houseBetaValues,
                                                              vnl_matrix <ScalarType> &qMatrix, unsigned int rank)
{
    unsigned int m = qrMatrix.rows();
    unsigned int n = qrMatrix.cols();

    if (m < n)
        return;

    qMatrix.set_size(m,m);
    qMatrix.set_identity();

    for (int j = rank - 1;j >= 0;--j)
    {
        unsigned int vecSize = m - j;
        vnl_vector <double> vtQ(vecSize);

        // Compute v^T B
        for (unsigned int i = 0;i < vecSize;++i)
        {
            vtQ[i] = qMatrix(j,i + j);
            for (unsigned int k = j + 1;k < m;++k)
                vtQ[i] += qrMatrix(k,j) * qMatrix(k,i + j);
        }

        for (unsigned int i = 0;i < vecSize;++i)
        {
            double constantValue = houseBetaValues[j];
            if (i != 0)
                constantValue *= qrMatrix(i + j,j);

            for (unsigned int k = 0;k < vecSize;++k)
                qMatrix(i + j,k + j) -= constantValue * vtQ[k];
        }
    }
}

} // end namespace anima
