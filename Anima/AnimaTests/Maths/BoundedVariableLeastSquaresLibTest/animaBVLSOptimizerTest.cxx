#include <animaBVLSOptimizer.h>
#include <itkTimeProbe.h>
#include <iostream>
#include <fstream>
#include <vnl/algo/vnl_qr.h>

int main()
{
    typedef anima::BVLSOptimizer OptimizerType;

    itk::TimeProbe tmpTime;
    tmpTime.Start();

    OptimizerType::Pointer optTest = OptimizerType::New();

    unsigned int dimX = 23;
    unsigned int dimY = 23;
    OptimizerType::MatrixType testData (dimX,dimY);
    std::ifstream dataMat("dataMatrix.txt");

    for (unsigned int i = 0;i < dimX;++i)
    {
        for (unsigned int j = 0;j < dimY;++j)
            dataMat >> testData(i,j);
    }

    std::cout << "Test data " << testData << std::endl;

    OptimizerType::ParametersType testPoints(dimX);
    std::ifstream bVector("bVector.txt");

    for (unsigned int i = 0;i < dimX;++i)
        bVector >> testPoints[i];

    std::cout << "Test points " << testPoints << std::endl;

    OptimizerType::ParametersType lowerBounds(dimY);
    std::ifstream lowBounds("lb.txt");

    for (unsigned int i = 0;i < dimY;++i)
        lowBounds >> lowerBounds[i];

    OptimizerType::ParametersType upperBounds(dimY);
    std::ifstream upBounds("ub.txt");

    for (unsigned int i = 0;i < dimY;++i)
        upBounds >> upperBounds[i];

    optTest->SetDataMatrix(testData);
    optTest->SetPoints(testPoints);

    optTest->SetLowerBounds(lowerBounds);
    optTest->SetUpperBounds(upperBounds);

    optTest->StartOptimization();

    tmpTime.Stop();

    std::cout << "Computation time: " << tmpTime.GetTotal() << std::endl;
    std::cout << "BVLS solution: " << optTest->GetCurrentPosition() << std::endl;
    std::cout << "BVLS residual: " << optTest->GetCurrentResidual() << std::endl;

    std::cout << "LS solution: " << vnl_qr <double> (testData).solve(testPoints) << std::endl;

    OptimizerType::ParametersType refOutput(dimY);
    std::ifstream refOutData("ref.txt");

    for (unsigned int i = 0;i < dimY;++i)
        refOutData >> refOutput[i];

//    std::cout << "Scipy solution: " << refOutput << std::endl;
//    std::cout << "Diff: " << refOutput - optTest->GetCurrentPosition() << std::endl;

    return EXIT_SUCCESS;
}
