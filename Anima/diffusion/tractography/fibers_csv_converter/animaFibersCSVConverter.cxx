#include <tclap/CmdLine.h>

#include <animaFibersWriter.h>
#include <animaFibersReader.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("Converts either any Anima-supported fiber file (VTK, VTP or FDS) to CSV or an input CSV file to Medinria FDS fiber file type. INRIA / IRISA - VisAGeS Team and Politecnico di Milano", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","input","input fiber file",true,"","input tracks",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }
    
    std::string inputFile = inArg.getValue();
    std::string inputExtension = inputFile.substr(inputFile.find_last_of('.') + 1);
    std::string inputName = inputFile.substr(0, inputFile.find_last_of('.'));
    
    std::string outputFile;
    if (inputExtension == "csv")
        outputFile = inputName + ".fds";
    else
        outputFile = inputName + ".csv";

    anima::FibersReader trackReader;
    trackReader.SetFileName(inputFile);
    trackReader.Update();

    anima::FibersWriter writer;
    writer.SetInputData(trackReader.GetOutput());
    writer.SetFileName(outputFile);
    std::cout << "Writing tracks: " << outputFile << std::endl;
    writer.Update();

    return EXIT_SUCCESS;
}
