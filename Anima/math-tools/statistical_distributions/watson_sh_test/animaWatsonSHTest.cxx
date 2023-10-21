#include <animaWatsonDistribution.h>
#include <iostream>

int main(int argc, char **argv)
{
    double kappa = std::stod(argv[1]);

    std::cout << "Testing SH approximation for kappa : " << kappa << std::endl;

    std::vector <double> coefs, derivatives;
    anima::WatsonDistribution watsonDistr;
    itk::Vector<double,3> meanAxis;
    meanAxis[0] = 0.0;
    meanAxis[1] = 0.0;
    meanAxis[2] = 1.0;
    watsonDistr.SetMeanAxis(meanAxis);
    watsonDistr.SetConcentrationParameter(kappa);
    watsonDistr.GetStandardWatsonSHCoefficients(coefs, derivatives);

    std::cout << "Dawson integral: " << anima::EvaluateDawsonIntegral(std::sqrt(kappa),false) << std::endl;

    std::cout << "Coefficients: ";
    for (unsigned int i = 0;i < coefs.size();++i)
        std::cout << coefs[i] / (4 * M_PI) << " ";
    std::cout << std::endl;

    std::cout << "Derivatives: ";
    for (unsigned int i = 0;i < coefs.size();++i)
        std::cout << derivatives[i] / (4 * M_PI) << " ";
    std::cout << std::endl;

    return EXIT_SUCCESS;
}
