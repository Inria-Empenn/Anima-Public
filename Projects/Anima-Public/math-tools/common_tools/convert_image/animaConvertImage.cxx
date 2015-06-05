#include <tclap/CmdLine.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkCommand.h>
#include <itkExtractImageFilter.h>

#include <animaReadWriteFunctions.h>
#include <animaRetrieveImageTypeMacros.h>

struct arguments
{
    std::string input, output;
    bool displayInfo;
};

template <class ImageType>
void
convert(const arguments &args)
{
    anima::writeImage<ImageType>(args.output, anima::readImage<ImageType>(args.input));
}

template <class ComponentType, int dimension>
void
checkIfComponentsAreVectors(itk::ImageIOBase::Pointer imageIO, const arguments &args)
{
    ANIMA_CHECK_IF_COMPONENTS_ARE_VECTORS(imageIO, ComponentType, dimension, convert, args)
}

template <class ComponentType>
void
retrieveNbDimensions(itk::ImageIOBase::Pointer imageIO, const arguments &args)
{
    ANIMA_RETRIEVE_NUMBER_OF_DIMENSIONS(imageIO, ComponentType, checkIfComponentsAreVectors, imageIO, args);
}

int main(int ac, const char** av)
{
    TCLAP::CmdLine cmd("animaConvertImage can be used to rewrite an image in a new compatible format ie: nii to nrrd. "
                       "It can also be used to display generic information on an image such as size, origin etc. using the -I option.");

    TCLAP::ValueArg<std::string> inputArg("i",
            "input",
            "Input image",
            true,
            "",
            "Input image.",
            cmd);

    TCLAP::ValueArg<std::string> outputArg("o",
            "output",
            "Output image",
            false,
            "",
            "Output image.",
            cmd);

    TCLAP::SwitchArg infoArg("I",
            "info",
            "Display info on the image or not",
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

    // Find out the type of the image in file
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inputArg.getValue().c_str(),
                                                                           itk::ImageIOFactory::ReadMode);

    if( !imageIO )
    {
        std::cerr << "Itk could not find suitable IO factory for the input" << std::endl;
        return EXIT_FAILURE;
    }

    // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
    imageIO->SetFileName(inputArg.getValue());
    imageIO->ReadImageInformation();

    arguments args;
    args.input = inputArg.getValue(); args.output = outputArg.getValue(); args.displayInfo = infoArg.getValue();

    int numDimensions = imageIO->GetNumberOfDimensions() - 1;
    if(args.displayInfo)
    {
        std::cout << "NUMBER OF DIMENSIONS: " << imageIO->GetNumberOfDimensions() << std::endl;

        std::cout << "SIZE : [";
        for(int dim = 0; dim < numDimensions; ++dim)
            std::cout << imageIO->GetDimensions(dim) << ", ";
        std::cout << imageIO->GetDimensions(numDimensions) << "]" << std::endl;

        std::cout << "ORIGIN : [";
        for(int dim = 0; dim < numDimensions; ++dim)
            std::cout << imageIO->GetOrigin(dim) << ", ";
        std::cout << imageIO->GetOrigin(numDimensions) << "]" << std::endl;

        std::cout << "SPACING : [";
        for(int dim = 0; dim < numDimensions; ++dim)
            std::cout << imageIO->GetSpacing(dim) << ", ";
        std::cout << imageIO->GetSpacing(numDimensions) << "]" << std::endl;

        std::cout << "DIRECTIONS :" << std::endl;
        for(unsigned int dir = 0; dir < imageIO->GetNumberOfDimensions(); ++dir)
        {
            std::cout << "\tDIM " << dir << " : [";
            for (int dim = 0; dim < numDimensions; ++dim)
                std::cout << imageIO->GetDirection(dim)[dir] << ", ";
            std::cout << imageIO->GetDirection(numDimensions)[dir] << "]" << std::endl;
        }
        std::cout << "NUMBER OF COMPONENTS: "<< imageIO->GetNumberOfComponents()<< std::endl;
        std::cout << "COMPONENT TYPE: "<< itk::ImageIOBase::GetComponentTypeAsString(imageIO->GetComponentType()) << std::endl;
    }

    if(args.output != "")
    {
        try
        {
            ANIMA_RETRIEVE_COMPONENT_TYPE(imageIO, retrieveNbDimensions, imageIO, args)
        }
        catch ( itk::ExceptionObject & err )
        {
            std::cerr << "Itk cannot convert, be sure to use valid arguments..." << std::endl;
            std::cerr << err << std::endl;
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}
