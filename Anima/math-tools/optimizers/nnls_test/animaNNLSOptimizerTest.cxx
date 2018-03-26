#include <animaNNLSOptimizer.h>
#include <itkTimeProbe.h>
#include <iostream>

int main()
{
    typedef anima::NNLSOptimizer OptimizerType;

    itk::TimeProbe tmpTime;
    tmpTime.Start();

    OptimizerType::Pointer optTest = OptimizerType::New();

    OptimizerType::MatrixType testData (3,3);
    testData(0,0) = 1;
    testData(0,1) = -3;
    testData(0,2) = 2;
    testData(1,0) = -3;
    testData(1,1) = 10;
    testData(1,2) = -5;
    testData(2,0) = 2;
    testData(2,1) = -5;
    testData(2,2) = 6;

    OptimizerType::ParametersType testPoints(3);
    testPoints[0] = 27;
    testPoints[1] = -78;
    testPoints[2] = 64;

    optTest->SetDataMatrix(testData);
    optTest->SetPoints(testPoints);

    optTest->StartOptimization();

    tmpTime.Stop();

    std::cout << "Computation time: " << tmpTime.GetTotal() << std::endl;
    std::cout << optTest->GetCurrentPosition() << std::endl;

    return 0;
}
