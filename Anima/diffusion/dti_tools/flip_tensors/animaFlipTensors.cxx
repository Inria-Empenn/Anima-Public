#include <animaFlipTensorsImageFilter.h>
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

    TCLAP::ValueArg<std::string> inArg("i", "input", "Input tensor image", true, "", "input tensor image", cmd);
    TCLAP::ValueArg<std::string> outArg("o", "output", "Output flipped tensor image", true, "", "output image", cmd);
    
    TCLAP::ValueArg<std::string> maskArg("m", "mask", "Computation mask", false, "", "mask", cmd);
    
    TCLAP::SwitchArg xDirArg("X", "xdir", "Flip X axis", cmd, false);
    TCLAP::SwitchArg yDirArg("Y", "ydir", "Flip Y axis", cmd, false);
    TCLAP::SwitchArg zDirArg("Z", "zdir", "Flip Z axis", cmd, false);

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

    typedef anima::FlipTensorsImageFilter<3> FilterType;
    typedef FilterType::InputImageType InputImageType;
    typedef FilterType::OutputImageType OutputImageType;
    typedef FilterType::MaskImageType MaskImageType;

    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(anima::readImage<InputImageType>(inArg.getValue()));
    filter->SetFlipXAxis(xDirArg.isSet());
    filter->SetFlipYAxis(yDirArg.isSet());
    filter->SetFlipZAxis(zDirArg.isSet());

    if (maskArg.getValue() != "")
        filter->SetComputationMask(anima::readImage<MaskImageType>(maskArg.getValue()));

    filter->SetNumberOfThreads(nbpArg.getValue());
    filter->AddObserver(itk::ProgressEvent(), callback );
    filter->Update();
    std::cout << std::endl;

    anima::writeImage<OutputImageType>(outArg.getValue(), filter->GetOutput());
    
    return EXIT_SUCCESS;
}// end of main
