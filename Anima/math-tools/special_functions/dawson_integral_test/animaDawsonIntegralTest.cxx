#include <animaErrorFunctions.h>

#include <itkTimeProbe.h>

#include <tclap/CmdLine.h>

#include <iostream>
#include <fstream>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ', ANIMA_VERSION);

    TCLAP::ValueArg<double> lbArg(
        "l", "lower-bound",
        "A numeric value specifying the lower bound of the domain (default: -5.0).",
        false, -5.0, "lower bound", cmd);
    TCLAP::ValueArg<double> ubArg(
        "u", "upper-bound",
        "A numeric value specifying the upper bound of the domain (default: +5.0).",
        false, 5.0, "upper bound", cmd);
    TCLAP::ValueArg<unsigned int> nPtsArg(
        "n", "nb-points",
        "An integer value specifying the size of the grid on which the domain is discretized (default: 1000).",
        false, 10000, "grid size", cmd);

    try
    {
        cmd.parse(argc, argv);
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    double lbValue = lbArg.getValue();
    double ubValue = ubArg.getValue();

    if (ubValue < lbValue)
    {
        std::cerr << "The lower bound should not exceed the upper bound." << std::endl;
        return EXIT_FAILURE;
    }

    unsigned int nbValues = nPtsArg.getValue();
    double stepValue = (ubValue - lbValue) / (static_cast<double>(nbValues) - 1.0);

    std::vector<double> inputValues(nbValues);
    std::vector<double> outputValues1(nbValues);
    std::vector<double> outputValues2(nbValues);
    std::vector<double> outputValues3(nbValues);

    double value = lbValue;
    for (unsigned int i = 0; i < nbValues; ++i)
    {
        inputValues[i] = value;
        value += stepValue;
    }

    itk::TimeProbe tmpTimer;

    tmpTimer.Start();

    try
    {
        for (unsigned int i = 0; i < nbValues; ++i)
            outputValues1[i] = anima::EvaluateDawsonIntegral(inputValues[i]);
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }

    tmpTimer.Stop();

    std::cout << "Dawson integral computed in " << tmpTimer.GetTotal() << " s" << std::endl;

    tmpTimer.Start();

    try
    {
        for (unsigned int i = 0; i < nbValues; ++i)
            outputValues2[i] = anima::EvaluateDawsonFunctionNR(inputValues[i]);
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }

    tmpTimer.Stop();

    std::cout << "EvaluateDawsonFunctionNR computed in " << tmpTimer.GetTotal() << " s" << std::endl;

    tmpTimer.Start();

    try
    {
        for (unsigned int i = 0; i < nbValues; ++i)
            outputValues3[i] = anima::EvaluateDawsonFunction(inputValues[i]);
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }

    tmpTimer.Stop();

    std::cout << "EvaluateDawsonFunction computed in " << tmpTimer.GetTotal() << " s" << std::endl;

    std::ofstream myfile;
    myfile.open("/Users/stamm-a/Downloads/example.csv");
    myfile << "x,y1,y2,y3\n";
    for (unsigned int i = 0; i < nbValues; ++i)
        myfile << std::to_string(inputValues[i]) << "," << std::to_string(outputValues1[i]) << "," << std::to_string(outputValues2[i]) << "," << std::to_string(outputValues3[i]) << "\n";
    myfile.close();

    return EXIT_SUCCESS;
}
