#include <animaFlipTensorImageFilter.h>
#include <tclap/CmdLine.h>
#include <animaReadWriteFunctions.h>
#include <animaRetrieveImageTypeMacros.h>

//Update progression of the process
void eventCallback(itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject *processObject = (itk::ProcessObject*) caller;
    std::cout << "\033[K\rProgression: " << (int)(processObject->GetProgress() * 100) << "%" << std::flush;
}

struct arguments
{
    std::string input, output, mask, axis;
    unsigned int nthreads;
};

template <class ComponentType, unsigned int ImageDimension>
void
flipTensors(itk::ImageIOBase::Pointer imageIO, const arguments &args)
{
    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);
    
    typedef anima::FlipTensorImageFilter<ComponentType,ImageDimension> FilterType;
    typedef typename FilterType::InputImageType InputImageType;
    typedef typename FilterType::OutputImageType OutputImageType;
    typedef typename FilterType::MaskImageType MaskImageType;
    
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput(anima::readImage<InputImageType>(args.input));
    filter->SetFlippedAxis(args.axis);
    
    if (args.mask != "")
        filter->SetComputationMask(anima::readImage<MaskImageType>(args.mask));
    
    filter->SetNumberOfThreads(args.nthreads);
    filter->AddObserver(itk::ProgressEvent(), callback);
    filter->Update();
    
    std::cout << std::endl;
    
    anima::writeImage<OutputImageType>(args.output, filter->GetOutput());
}

template <class ComponentType>
void
retrieveNbDimensions(itk::ImageIOBase::Pointer imageIO, const arguments &args)
{
    ANIMA_RETRIEVE_NUMBER_OF_DIMENSIONS(imageIO, ComponentType, flipTensors, imageIO, args)
}

int main(int argc, char **argv)
{

    TCLAP::CmdLine cmd("Flip tensors in a DTI volume.\nINRIA / IRISA - VisAGeS Team",' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i", "input", "Input tensor image.", true, "", "input image", cmd);
    TCLAP::ValueArg<std::string> outArg("o", "output", "Output tensor image.", true, "", "output image", cmd);
    
    TCLAP::ValueArg<std::string> maskArg("m", "mask", "Computation mask", false, "", "mask image", cmd);
    TCLAP::ValueArg<std::string> axisArg("a", "axis", "Axis to be flipped (choices are X, Y [default] or Z).", false, "Y", "axis name", cmd);

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
    
    // Retrieve image info
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inArg.getValue().c_str(),
                                                                           itk::ImageIOFactory::ReadMode);
    if (!imageIO)
    {
        std::cerr << "Itk could not find suitable IO factory for the input" << std::endl;
        return EXIT_FAILURE;
    }
    
    // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
    imageIO->SetFileName(inArg.getValue());
    imageIO->ReadImageInformation();
    
    arguments args;
    args.input = inArg.getValue();
    args.output = outArg.getValue();
    args.mask = maskArg.getValue();
    
    std::string axisStr = axisArg.getValue();
    std::transform(axisStr.begin(),axisStr.end(),axisStr.begin(),[](unsigned char c){ return std::tolower(c); });
    args.axis = axisStr;
    args.nthreads = nbpArg.getValue();
    
    try
    {
        ANIMA_RETRIEVE_COMPONENT_TYPE(imageIO, retrieveNbDimensions, imageIO, args);
    }
    catch (itk::ExceptionObject &err)
    {
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
    
}// end of main
