#include <animaExpTensorImageFilter.h>
#include <animaReadWriteFunctions.h>

#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team",' ', ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","inputlist","Log-tensors image",true,"","input log-tensor image",cmd);
    TCLAP::ValueArg<std::string> resArg("o","outputfile","Result tensor image",true,"","result tensor image",cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("T","nb-threads","Number of threads (default: all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"Number of threads",cmd);
    TCLAP::SwitchArg scaleArg("S","scale","The log-tensors have their non-diagonal terms scaled",cmd,false);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef anima::ExpTensorImageFilter <float> MainFilterType;

    MainFilterType::Pointer mainFilter = MainFilterType::New();

    mainFilter->SetNumberOfThreads(nbpArg.getValue());
    mainFilter->SetScaleNonDiagonal(scaleArg.isSet());
    mainFilter->SetInput(anima::readImage <MainFilterType::TInputImage> (inArg.getValue()));

    mainFilter->Update();

    std::cout << "Writing result to : " << resArg.getValue() << std::endl;

    anima::writeImage <MainFilterType::TOutputImage> (resArg.getValue(),mainFilter->GetOutput());

    return EXIT_SUCCESS;
}
