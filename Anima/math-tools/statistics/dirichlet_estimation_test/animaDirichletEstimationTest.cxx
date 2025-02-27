#include <animaDirichletDistribution.h>
#include <tclap/CmdLine.h>
#include <random>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ', ANIMA_VERSION);
    
    TCLAP::ValueArg<double> aArg("a", "a1", "First concentration parameter.", true, 1, "first concentration", cmd);
    TCLAP::ValueArg<double> bArg("b", "a2", "Second concentration parameter.", true, 1, "second concentration", cmd);
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
    
    std::vector<double> concentrationParameters(2);
    concentrationParameters[0] = aArg.getValue();
    concentrationParameters[1] = bArg.getValue();

    anima::DirichletDistribution dirichletDistribution;
    dirichletDistribution.SetConcentrationParameters(concentrationParameters);

    std::mt19937 generator(1234);

    std::vector <std::vector<double> > sample(nSampleArg.getValue());
    for (unsigned int i = 0;i < nSampleArg.getValue();++i)
        sample[i].resize(2);
    dirichletDistribution.Random(sample, generator);

    dirichletDistribution.Fit(sample, "mle");

    concentrationParameters = dirichletDistribution.GetConcentrationParameters();
    std::cout << "Concentration parameters: ";
    for (unsigned int i = 0;i < sample[0].size();++i)
        std::cout << concentrationParameters[i] << " ";
    std::cout << std::endl;

    return EXIT_SUCCESS;
}
