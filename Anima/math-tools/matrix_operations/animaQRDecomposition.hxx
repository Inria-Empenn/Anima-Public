#pragma once
#include "animaQRDecomposition.h"
#include <animaVectorOperations.h>
#include <limits>

namespace anima
{

template <typename ScalarType> void QRGivensDecomposition(vnl_matrix <ScalarType> &aMatrix, vnl_vector <ScalarType> &bVector)
{
    unsigned int m = aMatrix.rows();
    unsigned int n = aMatrix.cols();

    if (m < n)
        return;

    bool applyQToB = (bVector.size() == m);

    for (unsigned int j = 0;j < n;++j)
    {
        for (unsigned int i = m - 1;i >= j + 1;--i)
        {
            double bValue = aMatrix.get(i,j);
            double aValue = aMatrix.get(j,j);
            if (std::abs(bValue) <= std::numeric_limits <double>::epsilon())
                continue;

            // Compute Givens cos and sine values
            double cosValue, sinValue;
            double r = std::hypot(aValue,bValue);
            cosValue = aValue / r;
            sinValue = - bValue / r;

            for (unsigned int k = j;k < n;++k)
            {
                double ajkValue = aMatrix.get(j,k);
                double aikValue = aMatrix.get(i,k);
                double ajValue = cosValue * ajkValue - sinValue * aikValue;
                double aiValue = sinValue * ajkValue + cosValue * aikValue;

                aMatrix.put(j,k,ajValue);
                aMatrix.put(i,k,aiValue);
            }

            if (applyQToB)
            {
                double bjValue = cosValue * bVector[j] - sinValue * bVector[i];
                double biValue = sinValue * bVector[j] + cosValue * bVector[i];

                bVector[j] = bjValue;
                bVector[i] = biValue;
            }
        }
    }
}

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
        {
            double tmpVal = aMatrix.get(j,i);
            cVector[i] += tmpVal * tmpVal;
        }
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

    std::vector <double> housedVector;
    std::vector <double> housedTransposeA(n);
    double epsilon = std::numeric_limits<double>::epsilon();

    while (tau > 0.0)
    {
        ++rank;
        unsigned int r = rank - 1;
        unsigned int tmpIndex = pivotVector[k];
        pivotVector[k] = pivotVector[r];
        pivotVector[r] = tmpIndex;
        double tmp = cVector[k];
        cVector[k] = cVector[r];
        cVector[r] = tmp;

        if (r != k)
        {
            for (unsigned int i = 0;i < m;++i)
            {
                tmp = aMatrix(i,k);
                aMatrix.put(i,k,aMatrix.get(i,r));
                aMatrix.put(i,r,tmp);
            }
        }

        housedVector.resize(m - r);
        for (unsigned int i = r;i < m;++i)
            housedVector[i - r] = aMatrix.get(i,r);

        ComputeHouseholderVector(housedVector,houseBetaValues[r]);

        double diagonalValueTest = 0.0;
        for (unsigned int l = r;l < m;++l)
        {
            double workMatrixHouseValue = 0.0;
            unsigned int lIndex = l - r;
            if (lIndex != 0)
                workMatrixHouseValue = - houseBetaValues[r] * housedVector[0] * housedVector[lIndex];
            else
                workMatrixHouseValue = 1.0 - houseBetaValues[r] * housedVector[0] * housedVector[0];

            diagonalValueTest += workMatrixHouseValue * aMatrix.get(l,r);
        }

        if (r > 0)
        {
            if (std::abs(diagonalValueTest) < epsilon)
            {
                // cancel modifications so far and exit
                --rank;
                houseBetaValues[r] = 0.0;
                tmpIndex = pivotVector[k];
                pivotVector[k] = pivotVector[r];
                pivotVector[r] = tmpIndex;
                tmp = cVector[k];
                cVector[k] = cVector[r];
                cVector[r] = tmp;

                if (r != k)
                {
                    for (unsigned int i = 0;i < m;++i)
                    {
                        tmp = aMatrix.get(i,k);
                        aMatrix.put(i,k,aMatrix.get(i,r));
                        aMatrix.put(i,r,tmp);
                    }
                }

                tau = 0.0;
                continue;
            }
        }
        else
        {
            // Compute reference value for diagonal value test
            // taken from GSL rank test
            double basePower = std::floor(std::log(std::abs(diagonalValueTest)) / std::log(2.0));
            epsilon *= 20.0 * (m + n) * std::pow(2.0,basePower);
        }

        for (unsigned int i = 0;i < n - r;++i)
        {
            housedTransposeA[i] = 0.0;
            for (unsigned int j = 0;j < m - r;++j)
                housedTransposeA[i] += housedVector[j] * aMatrix.get(j + r, i + r);
        }

        for (unsigned int i = r;i < m;++i)
        {
            for (unsigned int j = r;j < n;++j)
            {
                double tmpValue = aMatrix.get(i,j) - houseBetaValues[r] * housedVector[i - r] * housedTransposeA[j - r];
                aMatrix.put(i,j,tmpValue);
            }
        }

        for (unsigned int i = rank;i < m;++i)
            aMatrix.put(i,r,housedVector[i - r]);

        for (unsigned int i = rank;i < n;++i)
        {
            double tmpVal = aMatrix.get(r,i);
            cVector[i] -= tmpVal * tmpVal;
        }

        tau = 0.0;

        if (rank < n)
        {
            k = rank;
            tau = cVector[rank];

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
}

template <typename ScalarType> void GetQtBFromQRPivotDecomposition(vnl_matrix <ScalarType> &qrMatrix, vnl_vector<ScalarType> &bVector,
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
            vtB += qrMatrix.get(i,j) * bVector[i];

        // bVector -= beta * v * vtB
        bVector[j] -= houseBetaValues[j] * vtB;

        for (unsigned int i = j + 1;i < m;++i)
            bVector[i] -= houseBetaValues[j] * qrMatrix.get(i,j) * vtB;
    }
}

template <typename ScalarType> void GetQMatrixQRPivotDecomposition(vnl_matrix <ScalarType> &qrMatrix, std::vector <ScalarType> &houseBetaValues,
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
                vtQ[i] += qrMatrix.get(k,j) * qMatrix.get(k,i + j);
        }

        for (unsigned int i = 0;i < vecSize;++i)
        {
            double constantValue = houseBetaValues[j];
            if (i != 0)
                constantValue *= qrMatrix.get(i + j,j);

            for (unsigned int k = 0;k < vecSize;++k)
            {
                double tmpVal = qMatrix.get(i + j,k + j) - constantValue * vtQ[k];
                qMatrix.put(i + j,k + j,tmpVal);
            }
        }
    }
}

} // end namespace anima
