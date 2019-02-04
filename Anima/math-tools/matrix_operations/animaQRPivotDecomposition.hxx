#pragma once
#include "animaQRPivotDecomposition.h"
#include <animaVectorOperations.h>

namespace anima
{

template <typename ScalarType> void QRPivotDecomposition(vnl_matrix <ScalarType> &aMatrix, std::vector <unsigned int> &transposePivotVector,
                                                         std::vector <ScalarType> &houseBetaValues, unsigned int &rank)
{
    unsigned int m = aMatrix.rows();
    unsigned int n = aMatrix.cols();

    if (m < n)
        return;

    transposePivotVector.resize(n);
    houseBetaValues.resize(n);

    for (unsigned int i = 0;i < n;++i)
    {
        transposePivotVector[i] = i;
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
        unsigned int tmpIndex = transposePivotVector[k];
        transposePivotVector[k] = transposePivotVector[r];
        transposePivotVector[r] = tmpIndex;
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

template <typename ScalarType> void GetQMatrixFromQRDecomposition(vnl_matrix <ScalarType> &qrMatrix, vnl_matrix <ScalarType> &qMatrix,
                                                                  std::vector <ScalarType> &houseBetaValues, unsigned int rank)
{
    unsigned int m = qrMatrix.rows();
    unsigned int n = qrMatrix.cols();

    if (m < n)
        return;

    qMatrix.set_size(m,m);
    qMatrix.set_identity();

    std::vector <ScalarType> workVector(m,0.0);
    for (int j = rank - 1;j >= 0;--j)
    {
        // Compute v^T Q
        for (unsigned int i = j;i < m;++i)
        {
            workVector[i] = qMatrix(j,i);
            for (unsigned int k = j + 1;k < m;++k)
                workVector[i] += qrMatrix(k,j) * qMatrix(k,i);
        }

        // Q -= beta * v * workVector
        for (unsigned int i = j;i < m;++i)
            qMatrix(j,i) -= houseBetaValues[j] * workVector[i];

        for (unsigned int i = j + 1;i < m;++i)
        {
            for (unsigned int k = j;k < m;++k)
                qMatrix(i,k) -= houseBetaValues[j] * qrMatrix(i,j) * workVector[k];
        }
    }
}

} // end namespace anima
