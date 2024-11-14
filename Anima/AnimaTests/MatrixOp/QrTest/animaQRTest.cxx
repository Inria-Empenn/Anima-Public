#include <animaQRDecomposition.h>
#include <iostream>
#include <fstream>

int main(int argc, char **argv)
{        
    unsigned int dimX = 12;
    unsigned int dimY = 6;

    vnl_matrix <double> testData (dimX,dimY);
    vnl_vector <double> bVector(dimX);
    bVector.fill(1.0);

    if (argc <= 1)
    {
        unsigned int pos = 10;
        for (unsigned int i = 0;i < dimX;++i)
        {
            for (unsigned int j = 0;j < dimY;++j)
                testData(i,j) = pos++;

            if (i != 0)
                testData(i,0) = 10;
        }
    }
    else
    {
        std::ifstream inputData(argv[1]);
        for (unsigned int i = 0;i < dimX;++i)
        {
            for (unsigned int j = 0;j < dimY;++j)
                inputData >> testData(i,j);
        }

        inputData.close();
    }

    std::cout << "Test data " << testData << std::endl;

    std::vector <double> betaValues(dimY,0.0);
    std::vector <unsigned int> pivotVector(dimY);
    unsigned int rank = 0;

    anima::QRPivotDecomposition(testData,pivotVector,betaValues,rank);

    std::cout.precision(20);
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

    anima::GetQtBFromQRPivotDecomposition(testData,bVector,betaValues,rank);

    std::cout << "QtB " << bVector << std::endl;

    vnl_matrix <double> qMatrix (dimX,dimX);
    anima::GetQMatrixQRPivotDecomposition(testData,betaValues,qMatrix,rank);

    std::cout << "Q matrix " << qMatrix << std::endl;

    return EXIT_SUCCESS;
}
