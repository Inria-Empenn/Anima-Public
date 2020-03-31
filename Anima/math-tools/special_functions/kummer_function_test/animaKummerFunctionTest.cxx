#include <iostream>
#include <tclap/CmdLine.h>
#include <animaKummerFunctions.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("Test Kummer function M. INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    TCLAP::ValueArg<double> xArg("x","input-value","Input value.",true,0,"input value",cmd);
    TCLAP::ValueArg<double> aArg("a","a-value","Input a value (default: -0.5).",false,-0.5,"input a value",cmd);
    TCLAP::ValueArg<double> bArg("b","b-value","Input b value (default:  1.0)",false,1,"input b value",cmd);
    
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    double xValue = xArg.getValue();
    double aValue = aArg.getValue();
    double bValue = bArg.getValue();

    std::cout << "M(" << xValue << ", ";
    std::cout << aValue << ", " << bValue << ") = ";
    std::cout << anima::GetKummerFunctionValue(xValue, aValue, bValue);
    std::cout << std::endl;

    std::cout << "exp(-" << xValue << ") M(" << xValue << ", ";
    std::cout << aValue << ", " << bValue << ") = ";
    std::cout << anima::GetScaledKummerFunctionValue(xValue, aValue, bValue);
    std::cout << std::endl;

    return EXIT_SUCCESS;
}
