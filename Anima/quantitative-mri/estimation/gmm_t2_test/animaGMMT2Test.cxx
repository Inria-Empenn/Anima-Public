#include <cmath>
#include <iostream>
#include <fstream>

#include <animaB1GMMRelaxometryCostFunction.h>

int main(int argc, char **argv)
{
    if (argc < 2)
        return EXIT_FAILURE;

    unsigned int numSignals = 32;
    anima::B1GMMRelaxometryCostFunction::ParametersType signals(numSignals);
    std::ifstream dataFile(argv[1]);
    for (unsigned int i = 0;i < numSignals;++i)
        dataFile >> signals[i];

    std::vector <double> means, variances;
    means.resize(3);
    means[0] = 20;
    means[1] = 100;
    means[2] = 2000;
    variances.resize(3);
    variances[0] = 25;
    variances[1] = 100;
    variances[2] = 6400;

    anima::B1GMMRelaxometryCostFunction::Pointer cost = anima::B1GMMRelaxometryCostFunction::New();
    cost->SetEchoSpacing(9);
    cost->SetExcitationFlipAngle(M_PI / 2.0);
    cost->SetT2RelaxometrySignals(signals);
    cost->SetGaussianMeans(means);
    cost->SetGaussianVariances(variances);
    cost->SetT1Value(1000);

    anima::B1GMMRelaxometryCostFunction::ParametersType b1(1);
    b1[0] = 1.1 * M_PI;
    std::cout << "Value for " << b1 << " " << cost->GetValue(b1) << std::endl;
    std::cout << cost->GetOptimalT2Weights() << std::endl;

    b1[0] = 1.2 * M_PI;
    std::cout << "Value for " << b1 << " " << cost->GetValue(b1) << std::endl;
    std::cout << cost->GetOptimalT2Weights() << std::endl;

    return EXIT_SUCCESS;
}
