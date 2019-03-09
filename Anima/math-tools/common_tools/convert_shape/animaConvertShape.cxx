#include <tclap/CmdLine.h>

#include <animaShapesReader.h>
#include <animaShapesWriter.h>

int main(int ac, const char** av)
{
    TCLAP::CmdLine cmd("animaConverthape can be used to convert a VTK, VTP, CSV or FDS (medInria fiber format) file into another format. It does not check before writing fds "
                       "that the file is actually made only of fibers.\n"
                       "INRIA / IRISA - VisAGeS/Empenn Team",' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inputArg("i","in","input data filename",true,"","input data",cmd);
    TCLAP::ValueArg<std::string> outputArg("o","out","output data filename",true,"","output data",cmd);

    try
    {
        cmd.parse(ac,av);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    anima::ShapesReader shapesReader;
    shapesReader.SetFileName(inputArg.getValue());
    try
    {
        shapesReader.Update();
    }
    catch(std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    anima::ShapesWriter shapesWriter;
    shapesWriter.SetInputData(shapesReader.GetOutput());
    shapesWriter.SetFileName(outputArg.getValue());

    try
    {
        shapesWriter.Update();
    }
    catch(std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
