#include <tclap/CmdLine.h>

#include <animaReadWriteFunctions.h>
#include <animaRetrieveImageTypeMacros.h>

#include <itkComposeImageFilter.h>

struct arguments
{
    std::vector<std::string> inputs;
    std::string output;
    bool fileList;
};

template <class InputImageType>
void
vectorize(const arguments &args)
{
    unsigned int nbComponents = args.inputs.size();
    typedef itk::VectorImage<typename InputImageType::PixelType, InputImageType::ImageDimension> OutputImageType;

    typedef itk::ComposeImageFilter<InputImageType, OutputImageType> ComposeImageFilterType;
    typename ComposeImageFilterType::Pointer composeImageFilter = ComposeImageFilterType::New();

    if (!args.fileList)
    {
        for(unsigned int i = 0; i < nbComponents; ++i)
            composeImageFilter->SetInput(i, anima::readImage<InputImageType>(args.inputs[i]));
    }
    else
    {
        std::string inputName = args.inputs[0];
        anima::setMultipleImageFilterInputsFromFileName<InputImageType,ComposeImageFilterType>(inputName,composeImageFilter);
    }

    composeImageFilter->Update();

    anima::writeImage<OutputImageType>(args.output, composeImageFilter->GetOutput());
}

template <class ComponentType, int dimension>
void
checkIfComponentsAreVectors(itk::ImageIOBase::Pointer imageIO, const arguments &args)
{
    if(imageIO->GetNumberOfComponents() != 1)
    {
        itk::ExceptionObject excp(__FILE__, __LINE__, "Number of components not supported.", ITK_LOCATION);
        throw excp;
    }
    else
    {
        typedef itk::Image<ComponentType, dimension> ImageType;
        vectorize<ImageType>(args);
    }
}

template <class ComponentType>
void
retrieveNbDimensions(itk::ImageIOBase::Pointer imageIO, const arguments &args)
{
    ANIMA_RETRIEVE_NUMBER_OF_DIMENSIONS(imageIO, ComponentType, checkIfComponentsAreVectors, imageIO, args);
}

int main(int ac, const char** av)
{
    TCLAP::CmdLine cmd("animaVectorizeImage is used to combine several scalar images into a multicomponent Vector image."
                       "INRIA / IRISA - VisAGeS/Empenn Team",
                       ' ',
                       ANIMA_VERSION);

    TCLAP::MultiArg<std::string> inputArg("i",
            "inputs",
            "input images (list in text file or multiple arguments)",
            true,
            "input images",
            cmd);

    TCLAP::ValueArg<std::string> outputArg("o",
            "output",
            "output image",
            true,
            "",
            "output image",
            cmd);

    try
    {
        cmd.parse(ac,av);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inputArg.getValue()[0].c_str(),itk::ImageIOFactory::ReadMode);

    arguments args;
    args.fileList = false;

    if( !imageIO )
    {
        std::ifstream fileIn (inputArg.getValue()[0].c_str());
        char tmpStr[2048];
        fileIn.getline(tmpStr,2048);

        while ((strcmp(tmpStr,"") == 0) && (!fileIn.eof()))
            fileIn.getline(tmpStr,2048);

        if (!fileIn.eof())
            imageIO = itk::ImageIOFactory::CreateImageIO(tmpStr,itk::ImageIOFactory::ReadMode);

        fileIn.close();

        if( !imageIO )
        {
            std::cerr << "Itk could not find suitable IO factory for the input" << std::endl;
            return EXIT_FAILURE;
        }

        imageIO->SetFileName(tmpStr);
        imageIO->ReadImageInformation();
        args.fileList = true;
    }

    // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
    if (!args.fileList)
    {
        imageIO->SetFileName(inputArg.getValue()[0]);
        imageIO->ReadImageInformation();
    }

    args.inputs = inputArg.getValue(); args.output = outputArg.getValue();

    if(args.output != "")
    {
        try
        {
            ANIMA_RETRIEVE_COMPONENT_TYPE(imageIO, retrieveNbDimensions, imageIO, args)
        }
        catch ( itk::ExceptionObject & err )
        {
            std::cerr << "Itk cannot vectorize, be sure to use valid arguments..." << std::endl;
            std::cerr << err << std::endl;
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}
