#include <iostream>
#include <tclap/CmdLine.h>
#include <animaKummerFunctions.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("Test Kummer function M. INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    TCLAP::ValueArg<double> inArg("i","input-value","Input value",true,0,"input value",cmd);
    
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }
    
    std::cout << "Input value: " << inArg.getValue() << std::endl;
    std::cout << "Kummer M value via direct ARB: " << anima::KummerFunction(inArg.getValue(), -0.5, 1.0) << std::endl;
    std::cout << "Kummer M value via exponentially-scaled modified Bessel functions: " << anima::OneHalfLaguerreFunction(inArg.getValue()) << std::endl;

    return EXIT_SUCCESS;
}
