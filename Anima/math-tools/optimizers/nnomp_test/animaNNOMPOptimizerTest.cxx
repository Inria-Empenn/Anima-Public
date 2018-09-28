#include <animaNNOrthogonalMatchingPursuitOptimizer.h>
#include <itkTimeProbe.h>
#include <iostream>
#include <fstream>

int main()
{
    typedef anima::NNOrthogonalMatchingPursuitOptimizer OptimizerType;

    itk::TimeProbe tmpTime;
    tmpTime.Start();

    OptimizerType::Pointer optTest = OptimizerType::New();

    OptimizerType::MatrixType testData (288,362);
    std::ifstream inputDataMatrix("dataMatrix.txt");

    for (unsigned int i = 0;i < 288;++i)
    {
        for (unsigned int j = 0;j < 362;++j)
            inputDataMatrix >> testData(i,j);
    }

    inputDataMatrix.close();
    OptimizerType::ParametersType testPoints(288);
    std::ifstream inputDataSignals("signals.txt");

    for (unsigned int i = 0;i < 288;++i)
        inputDataSignals >> testPoints[i];

    inputDataSignals.close();

    optTest->SetDataMatrix(testData);
    optTest->SetPoints(testPoints);

    optTest->SetMaximalNumberOfWeights(3);
    optTest->SetIgnoredIndexesUpperBound(2);

    optTest->StartOptimization();

    tmpTime.Stop();

    std::cout << "Computation time: " << tmpTime.GetTotal() << std::endl;
    std::cout << optTest->GetCurrentPosition() << std::endl;

    return EXIT_SUCCESS;
}
