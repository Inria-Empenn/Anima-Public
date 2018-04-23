#include <tclap/CmdLine.h>

#include <animaShapesReader.h>
#include <animaShapesWriter.h>

int main(int ac, const char** av)
{
    TCLAP::CmdLine cmd("animaConverthapes can be used to rewrite a VTK, CSV or FDS (fiber format) file into another format. It does not check though before writing fds that the file is "
                       "actually made only of fibers.\n"
                       "INRIA / IRISA - VisAGeS Team",' ',ANIMA_VERSION);

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
