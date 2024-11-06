#include <animaGammaDistribution.h>
#include <tclap/CmdLine.h>
#include <random>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ', ANIMA_VERSION);
    
    TCLAP::ValueArg<double> kappaArg("k", "kappa", "Shape parameter.", true, 1, "shape parameter", cmd);
    TCLAP::ValueArg<double> thetaArg("t", "theta", "Scale parameter.", true, 1, "scale parameter", cmd);
    TCLAP::ValueArg<unsigned int> nSampleArg("n", "nsample", "Sample size.", false, 10000, "sample size", cmd);
    TCLAP::ValueArg<std::string> methodArg("", "method", "Estimation method. Choices are: mle, biased-closed-form or unbiased-closed-form", false, "mle", "estimation method", cmd);
    
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }
    
    double kappa = kappaArg.getValue();
    double theta = thetaArg.getValue();

    anima::GammaDistribution gammaDistribution;
    gammaDistribution.SetShapeParameter(kappa);
    gammaDistribution.SetScaleParameter(theta);

    std::mt19937 generator(1234);

    std::vector<double> sample(nSampleArg.getValue());
    gammaDistribution.Random(sample, generator);

    gammaDistribution.Fit(sample, methodArg.getValue());

    std::cout << "Shape parameter: " << gammaDistribution.GetShapeParameter() << std::endl;
    std::cout << "Scale parameter: " << gammaDistribution.GetScaleParameter() << std::endl;

    return EXIT_SUCCESS;
}
