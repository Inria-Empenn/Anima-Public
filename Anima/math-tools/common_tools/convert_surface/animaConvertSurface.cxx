#include <tclap/CmdLine.h>

#include <animaFibersReader.h>
#include <animaFibersWriter.h>

int main(int ac, const char** av)
{
    TCLAP::CmdLine cmd("animaConvertSurface can be used to rewrite a VTK or FDS (fiber format) file into another format. It does not check though before writing fds that the file is "
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

    anima::FibersReader fibersReader;
    fibersReader.SetFileName(inputArg.getValue());
    try
    {
        fibersReader.Update();
    }
    catch(std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    anima::FibersWriter fibersWriter;
    fibersWriter.SetInputData(fibersReader.GetOutput());
    fibersWriter.SetFileName(outputArg.getValue());

    try
    {
        fibersWriter.Update();
    }
    catch(std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
