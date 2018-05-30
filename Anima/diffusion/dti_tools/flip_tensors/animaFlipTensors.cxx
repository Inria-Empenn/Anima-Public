#include <animaFlipTensorImageFilter.h>
#include <tclap/CmdLine.h>
#include <animaReadWriteFunctions.h>

//Update progression of the process
void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
    std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
}

int main(int argc, char **argv)
{

    TCLAP::CmdLine cmd("Flip tensors in a DTI volume.\nINRIA / IRISA - VisAGeS Team",' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i", "input", "Input tensor image.", true, "", "input image", cmd);
    TCLAP::ValueArg<std::string> outArg("o", "output", "Output tensor image.", true, "", "output image", cmd);
    
    TCLAP::ValueArg<std::string> maskArg("m", "mask", "Computation mask", false, "", "mask image", cmd);
    TCLAP::ValueArg<std::string> axisArg("a", "axis", "Axis to be flipped (Choices are none [default], X, Y or Z).", false, "", "axis name", cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("p", "nthreads", "Number of thread to use (default: all)", false, itk::MultiThreader::GetGlobalDefaultNumberOfThreads(), "number of thread", cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);

    typedef anima::FlipTensorImageFilter<3> FilterType;
    typedef FilterType::InputImageType InputImageType;
    typedef FilterType::OutputImageType OutputImageType;
    typedef FilterType::MaskImageType MaskImageType;

    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(anima::readImage<InputImageType>(inArg.getValue()));
    filter->SetFlippedAxis(axisArg.getValue());

    if (maskArg.getValue() != "")
        filter->SetComputationMask(anima::readImage<MaskImageType>(maskArg.getValue()));

    filter->SetNumberOfThreads(nbpArg.getValue());
    filter->AddObserver(itk::ProgressEvent(), callback );
    
    try
    {
        filter->Update();
    }
    catch (itk::ExceptionObject &err)
    {
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }
    
    std::cout << std::endl;

    anima::writeImage<OutputImageType>(outArg.getValue(), filter->GetOutput());
    
    return EXIT_SUCCESS;
    
}// end of main
