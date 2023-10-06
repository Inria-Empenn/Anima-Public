#include <cmath>

#include <animaDDIDistribution.h>
#include <animaHyperbolicFunctions.h>
#include <vnl/vnl_matrix.h>

#include <tclap/CmdLine.h>

int main(int argc,  char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<double> RadiusArg("r","radius","Radius",false,1,"Radius",cmd);
    TCLAP::ValueArg<double> XDirArg("x","xdir","X coord of orientation",false,1,"X coord of orientation",cmd);
    TCLAP::ValueArg<double> YDirArg("y","ydir","Y coord of orientation",false,0,"Y coord of orientation",cmd);
    TCLAP::ValueArg<double> ZDirArg("z","zdir","Z coord of orientation",false,0,"Z coord of orientation",cmd);
    TCLAP::ValueArg<double> KappaArg("k","kappa","Concentration",false,10,"Concentration",cmd);
    TCLAP::ValueArg<double> MeanDiffArg("d","mean-diff","Mean diffusivity",false,0,"Mean diffusivity",cmd);
    TCLAP::ValueArg<double> IntraVArg("i","intra-vol","Intra-axonal volume fraction",false,0,"Intra-axonal volume fraction",cmd);
       
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }
    
    double step = 1e-2;
    double inStep = 1e-2;
    unsigned int gridSize = (unsigned int)(1.0/step+1);
    std::vector <double> theta(gridSize,0);
    std::vector <double> phi(gridSize,0);
    std::vector <double> inVec(3,0);
    std::vector <double> muVec(3,0);
    muVec[0] = XDirArg.getValue();
    muVec[1] = YDirArg.getValue();
    muVec[2] = ZDirArg.getValue();
    std::vector <double> integrand((unsigned int)(1.0/inStep+1),0);
    vnl_matrix <double> prob;
    prob.set_size(gridSize,gridSize);
    prob.fill(0);
    
    double k = KappaArg.getValue();
    double d = MeanDiffArg.getValue();
    double n = IntraVArg.getValue();
    double dFiber = 1.71e-3 / (1.0 - 2.0 * n * anima::xi(k));
    if (d == 0)
        d = dFiber;
    double r = sqrt(RadiusArg.getValue() * dFiber);
    
    for (unsigned int i = 0;i < gridSize;++i)
    {
        theta[i] = i*step*M_PI;
        inVec[2] = r * cos(theta[i]);
        for (unsigned int j = 0;j < gridSize;++j)
        {
            phi[j] = j*step*2.0*M_PI;
            inVec[0] = r * sin(theta[i]) * cos(phi[j]);
            inVec[1] = r * sin(theta[i]) * sin(phi[j]);
            prob(i,j) = anima::ComputeSymmetricPDF(inVec, muVec, k, d, n, inStep, integrand);
        }
    }
    
    std::cout << "List of Theta values: " << std::endl;
    for (unsigned int i = 0;i < gridSize;++i)
        std::cout << theta[i] << " ";
    std::cout << std::endl;

    std::cout << "List of Phi values: " << std::endl;
    for (unsigned int i = 0;i < gridSize;++i)
        std::cout << phi[i] << " ";
    std::cout << std::endl;
    
    std::cout << "List of Proba values: " << std::endl;
    std::cout << prob << std::endl;
    
    return 0;
}
