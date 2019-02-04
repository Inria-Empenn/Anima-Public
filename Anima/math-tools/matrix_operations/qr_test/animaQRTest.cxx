#include <animaQRPivotDecomposition.h>
#include <iostream>

int main()
{
    unsigned int dimX = 4;
    unsigned int dimY = 3;
    vnl_matrix <double> testData (dimX,dimY);

    unsigned int pos = 1;
    for (unsigned int i = 0;i < dimX;++i)
    {
        for (unsigned int j = 0;j < dimY;++j)
            testData(i,j) = pos++;

        if (i != 1)
            testData(i,0) = 1;
    }

    std::cout << "Test data " << testData << std::endl;

    std::vector <double> betaValues(dimY,0.0);
    std::vector <unsigned int> pivotVector(dimY);
    unsigned int rank = 0;

    anima::QRPivotDecomposition(testData,pivotVector,betaValues,rank);

    std::cout << "Matrix rank " << rank << std::endl;
    std::cout << "Transformed matrix " << testData << std::endl;

    std::cout << "Pivot vector ";
    for (unsigned int i = 0;i < pivotVector.size();++i)
        std::cout << pivotVector[i] << " ";
    std::cout << std::endl;

    std::cout << "Beta Values ";
    for (unsigned int i = 0;i < betaValues.size();++i)
        std::cout << betaValues[i] << " ";
    std::cout << std::endl;

    vnl_matrix <double> qMatrix;
    anima::GetQMatrixFromQRDecomposition(testData,qMatrix,betaValues,rank);

    std::cout << "Q matrix " << qMatrix << std::endl;

    return EXIT_SUCCESS;
}
