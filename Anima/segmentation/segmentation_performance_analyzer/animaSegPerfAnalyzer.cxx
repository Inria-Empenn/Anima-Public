#include "SegPerfApp.h"
#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    CSegPerfApp oSegPerfApp;

    bool needWork = false;
    try
    {
        needWork = oSegPerfApp.init(argc, argv);
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    if (!oSegPerfApp.checkParamsCoherence())
        return EXIT_FAILURE;

    if (!needWork)
        return EXIT_SUCCESS;

    oSegPerfApp.checkOutputCoherence();
    oSegPerfApp.prepareOutput();

    try
    {
        oSegPerfApp.play();
    }
    catch(std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
