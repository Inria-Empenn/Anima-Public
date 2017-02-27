#include <animaGeneralizedFAImageFilter.h>
#include <animaReadWriteFunctions.h>
#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","inputodf","ODF volume",true,"","ODF volume",cmd);
    TCLAP::ValueArg<std::string> resArg("o","output","Result image",true,"","result GFA image",cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default: all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef anima::GeneralizedFAImageFilter <float> MainFilterType;

    MainFilterType::Pointer mainFilter = MainFilterType::New();
    mainFilter->SetInput(anima::readImage <MainFilterType::TInputImage> (inArg.getValue()));
    mainFilter->SetNumberOfThreads(nbpArg.getValue());

    mainFilter->Update();

    anima::writeImage <MainFilterType::TOutputImage> (resArg.getValue(),mainFilter->GetOutput());

    return EXIT_SUCCESS;
}
