#include <animaNNLSOptimizer.h>
#include <itkTimeProbe.h>
#include <iostream>
#include <fstream>

int main()
{
    typedef anima::NNLSOptimizer OptimizerType;

    itk::TimeProbe tmpTime;
    tmpTime.Start();

    OptimizerType::Pointer optTest = OptimizerType::New();

    unsigned int dimX = 72;
    unsigned int dimY = 2;
    OptimizerType::MatrixType testData (dimX,dimY);
    std::ifstream dataMat("/Users/ocommowi/Documents/Tmp2/dataMatrix.txt");

    for (unsigned int i = 0;i < dimX;++i)
    {
        for (unsigned int j = 0;j < dimY;++j)
            dataMat >> testData(i,j);
    }

    std::cout << "Test data " << testData << std::endl;

    OptimizerType::ParametersType testPoints(dimX);
    std::ifstream bVector("/Users/ocommowi/Documents/Tmp2/bVector.txt");

    for (unsigned int i = 0;i < dimX;++i)
        bVector >> testPoints[i];

    std::cout << "Test points " << testPoints << std::endl;

    optTest->SetDataMatrix(testData);
    optTest->SetPoints(testPoints);
    optTest->SetSquaredProblem(false);

    optTest->StartOptimization();

    tmpTime.Stop();

    std::cout << "Computation time: " << tmpTime.GetTotal() << std::endl;
    std::cout << optTest->GetCurrentPosition() << std::endl;

    return EXIT_SUCCESS;
}
