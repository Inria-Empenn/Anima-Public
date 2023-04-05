#include <animaDirichletDistribution.h>
#include <tclap/CmdLine.h>
#include <random>
#include <vnl/vnl_matrix.h>

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

    vnl_matrix<double> sample(nSampleArg.getValue(), 2);
    dirichletDistribution.Random(sample, generator);

    std::cout << sample << std::endl;

    for (unsigned int i = 0;i < sample.rows();++i)
    {
        double sumValue = 0.0;
        for (unsigned int j = 0;j < sample.cols();++j)
            sumValue += sample(i, j);
        std::cout << sumValue << " ";
    }
    std::cout << std::endl;

    return EXIT_SUCCESS;
}
