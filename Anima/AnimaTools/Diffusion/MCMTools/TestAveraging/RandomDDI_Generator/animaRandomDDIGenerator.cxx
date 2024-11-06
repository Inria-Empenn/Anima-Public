#include <cmath>
#include <fstream>
#include <iostream>

#include<cstdlib>
#include<ctime>
#include <animaMatrixOperations.h>

#include <tclap/CmdLine.h>

using namespace anima;

int main(int argc, char *argv[])
{

    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> firstNameArg("f","firstName","name of the output file with .txt",true,"","name of the output file with the name of 4 files",cmd);
    TCLAP::ValueArg<std::string> pixelNameArg("s","secondName","name of the output files without .txt",true,"","name of the  4 output files",cmd);

    TCLAP::ValueArg<double> kappaLow("","kappaLow","kappa lwo bound",false,0,"kappa lwo bound",cmd);
    TCLAP::ValueArg<double> kappaHigh("","kappaHigh","kappa high bound",false,20,"kappa high bound",cmd);

    TCLAP::ValueArg<double> vLow("","vLow","v lwo bound",false,0,"v lwo bound",cmd);
    TCLAP::ValueArg<double> vHigh("","vHigh","v high bound",false,1,"v high bound",cmd);

    TCLAP::ValueArg<double> dLow("","dLow","d lwo bound",false,0.0005,"d lwo bound",cmd);
    TCLAP::ValueArg<double> dHigh("","dHigh","d high bound",false,0.005,"d high bound",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    double kappa, v, d, theta, phi;
    vnl_vector<double> angles(3, 0), direction(3, 0);

    srand(time(NULL));

    for (int i=0; i<4; i++)
    {

        kappa = rand()/(double)RAND_MAX * (kappaHigh.getValue() - kappaLow.getValue()) + kappaLow.getValue();
        v = rand()/(double)RAND_MAX * (vHigh.getValue() - vLow.getValue()) + vLow.getValue();
        d = rand()/(double)RAND_MAX * (dHigh.getValue() - dLow.getValue()) + dLow.getValue();
        theta = rand()/(double)RAND_MAX  * 2 * M_PI;
        phi =  rand()/(double)RAND_MAX  * M_PI;

        angles(0) = phi;
        angles(1) = theta;
        angles(2) = 1;
        TransformSphericalToCartesianCoordinates(angles, direction);

        std::stringstream convertNumber;
        convertNumber << i;

        std::string namePixel = pixelNameArg.getValue() + convertNumber.str() + ".txt";
        std::ofstream fichierPixel;
        fichierPixel.open(namePixel.c_str());

        if (fichierPixel)
        {
            fichierPixel << direction(0) << std::endl;
            fichierPixel << direction(1) << std::endl;
            fichierPixel << direction(2) << std::endl;
            fichierPixel << kappa << std::endl;
            fichierPixel << d << std::endl;
            fichierPixel << v << std::endl;
            fichierPixel << 1.0 << std::endl;

            fichierPixel.close();
        }
        else
            std::cerr << "Erreur à l'ouverture !" << std::endl;
    }

    std::ofstream fichier;
    fichier.open(firstNameArg.getValue().c_str());

    for (int i=0; i<4; i++)
        fichier << pixelNameArg.getValue() << i << ".txt" << std::endl;

    fichier.close();

    return 0;
}
