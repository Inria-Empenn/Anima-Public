#include <animaWatsonDistribution.h>
#include <tclap/CmdLine.h>
#include <random>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ', ANIMA_VERSION);
    
    TCLAP::ValueArg<double> kappaArg("k", "kappa", "Shape parameter.", true, 1, "shape parameter", cmd);
    TCLAP::ValueArg<unsigned int> nSampleArg("n", "nsample", "Sample size.", false, 10000, "sample size", cmd);
    
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
    
    anima::WatsonDistribution watsonDistribution;

    using VectorType = itk::Vector<double,3>;
    VectorType meanAxis;
    meanAxis.Fill(1.0 / std::sqrt(3.0));
    double normValue = meanAxis.GetNorm();
    meanAxis /= normValue;
    watsonDistribution.SetMeanAxis(meanAxis);
    watsonDistribution.SetConcentrationParameter(kappa);

    std::mt19937 generator(1234);

    std::vector<VectorType> sample(nSampleArg.getValue());
    watsonDistribution.Random(sample, generator);

    watsonDistribution.Fit(sample, "");

    std::cout << "----- Mean Axis -----" << std::endl;
    std::cout << "- Expected : " << meanAxis << std::endl;
    std::cout << "- Estimated: " << watsonDistribution.GetMeanAxis() << std::endl;

    std::cout << std::endl;

    std::cout << "--- Concentration ---" << std::endl;
    std::cout << "- Expected : " << kappa << std::endl;
    std::cout << "- Estimated: " << watsonDistribution.GetConcentrationParameter() << std::endl;
    
    return EXIT_SUCCESS;
}
