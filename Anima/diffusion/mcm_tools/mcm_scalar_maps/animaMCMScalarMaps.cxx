#include <animaMCMFileReader.h>
#include <animaReadWriteFunctions.h>
#include <animaMCMScalarMapsImageFilter.h>

#include <itkTimeProbe.h>
#include <tclap/CmdLine.h>

int main(int argc,  char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","ddi","DDI volume",true,"","DDI volume",cmd);
    TCLAP::ValueArg<std::string> outFAArg("f","out-fa","Result FA image",true,"","result FA image",cmd);
    TCLAP::ValueArg<std::string> outMDArg("m","out-md","Result MD image",true,"","result MD image",cmd);

    TCLAP::ValueArg<std::string> outIsoRWArg("r","out-iso-r-w","Result iso restricted weight image",false,"","result iso R image",cmd);
    TCLAP::ValueArg<std::string> outAnisoRWArg("a","out-aniso-r-w","Result anisotropic weight image",false,"","result anisotropic image",cmd);
    TCLAP::ValueArg<std::string> outFWArg("F","out-fw","Result free water weight image",false,"","result FW image",cmd);

    TCLAP::SwitchArg includeIsoArg("I", "inc-iso", "Include isotropic contributions in FA and MD?", cmd, false);

    TCLAP::ValueArg<unsigned int> nbThreadsArg("T", "nb-threads", "Number of threads to run on (default: all cores)", false, itk::MultiThreader::GetGlobalDefaultNumberOfThreads(), "number of threads", cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef double ScalarType;
    typedef anima::MCMScalarMapsImageFilter <ScalarType> MainFilterType;

    MainFilterType::Pointer mainFilter = MainFilterType::New();

    anima::MCMFileReader <double, 3> mcmReader;
    mcmReader.SetFileName(inArg.getValue());
    mcmReader.Update();

    mainFilter->SetInput(mcmReader.GetModelVectorImage());
    mainFilter->SetIncludeIsotropicWeights(includeIsoArg.isSet());
    mainFilter->SetNumberOfThreads(nbThreadsArg.getValue());
    mainFilter->Update();

    anima::writeImage <MainFilterType::OutputImageType> (outFAArg.getValue(),mainFilter->GetOutput(3));
    anima::writeImage <MainFilterType::OutputImageType> (outMDArg.getValue(),mainFilter->GetOutput(4));

    if (outFWArg.getValue() != "")
        anima::writeImage <MainFilterType::OutputImageType> (outFWArg.getValue(),mainFilter->GetOutput(0));

    if (outIsoRWArg.getValue() != "")
        anima::writeImage <MainFilterType::OutputImageType> (outIsoRWArg.getValue(),mainFilter->GetOutput(1));

    if (outAnisoRWArg.getValue() != "")
        anima::writeImage <MainFilterType::OutputImageType> (outAnisoRWArg.getValue(),mainFilter->GetOutput(2));

    return EXIT_SUCCESS;
}
