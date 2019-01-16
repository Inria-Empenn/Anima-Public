#pragma once
#include "animaQRPivotDecomposition.h"
#include <animaVectorOperations.h>

namespace anima
{

template <typename ScalarType> void QRPivotDecomposition(vnl_matrix <ScalarType> &aMatrix, std::vector <unsigned int> &pivotVector, unsigned int &rank)
{
    unsigned int m = aMatrix.rows();
    unsigned int n = aMatrix.cols();

    if (m < n)
        return;

    pivotVector.resize(n);
    for (unsigned int i = 0;i < n;++i)
        pivotVector[i] = i;

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
    double beta = 0.0;
    vnl_matrix <double> workMatrixHouse(m,m);
    vnl_matrix <double> workMatrix(m,n);

    while (tau > 0)
    {
        ++rank;
        unsigned int r = rank - 1;
        pivotVector[r] = k;
        double tmp = cVector[k];
        cVector[k] = cVector[r];
        cVector[r] = tmp;

        for (unsigned int i = 0;i < m;++i)
        {
            tmp = aMatrix(i,k);
            aMatrix(i,k) = aMatrix(i,r);
            aMatrix(i,r) = tmp;
        }

        housedVector.resize(m - r + 1);
        for (unsigned int i = r;i < m;++i)
            housedVector[i - r] = aMatrix(i,r);

        ComputeHouseholderVector(housedVector,vVector,beta);

        for (unsigned int i = r;i < m;++i)
        {
            for (unsigned int j = r;j < n;++j)
                workMatrix(i,j) = aMatrix(i,j);
        }

        for (unsigned int i = 0;i < m - r + 1;++i)
        {
            for (unsigned int j = i;j < m - r + 1;++j)
            {
                workMatrixHouse(i,j) = - beta * vVector[i] * vVector[j];
                if (i != j)
                    workMatrixHouse(j,i) = workMatrixHouse(i,j);
                else
                    workMatrixHouse(i,j) += 1.0;
            }
        }

        for (unsigned int i = 0;i < m - r + 1;++i)
        {
            unsigned int iIndex = i + r;
            for (unsigned int j = r;j < n;++j)
            {
                aMatrix(iIndex,j) = 0.0;
                for (unsigned int l = 0;l < m - r + 1;++l)
                {
                    unsigned int lIndex = l + r;
                    aMatrix(iIndex,j) += workMatrixHouse(i,l) * workMatrix(lIndex,j);
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
        }

        for (unsigned int i = rank + 1;i < n;++i)
        {
            if (tau < cVector[i])
            {
                k = i;
                tau = cVector[i];
            }
        }
    }
}

} // end namespace anima
