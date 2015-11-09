#include <tclap/CmdLine.h>

#include <itkImage.h>
#include <itkOrientImageFilter.h>

#include <animaReorientation.h>
#include <animaReadWriteFunctions.h>
#include <animaRetrieveImageTypeMacros.h>

struct arguments
{
    std::string input, output, reorient, space;
    bool displayInfo;
};

template <class ImageType>
void
convert(const arguments &args)
{
    typename ImageType::Pointer input = anima::readImage<ImageType>(args.input);

    if(args.space != "")
    {
        typename ImageType::Pointer spaceImage = anima::readImage<ImageType>(args.space);

        input->CopyInformation(spaceImage);
    }

    if(args.reorient == "AXIAL")
        input = anima::reorientImage<ImageType>(input,
                                                itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
    else if(args.reorient == "CORONAL")
        input = anima::reorientImage<ImageType>(input,
                                                itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSA);
    else if(args.reorient == "SAGITTAL")
        input = anima::reorientImage<ImageType>(input,
                                                itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASL);

    input->SetMetaDataDictionary(itk::MetaDataDictionary());
    anima::writeImage<ImageType>(args.output, input);
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
                       "It can also be used to display generic information on an image such as size, origin etc. using the --info option."
                       "If an image file is given to the --space option it will be used to rewrite the input in the same"
                       "real coordinates repair."
                       "The --reorient option allow to reorient the image in either the AXIAL, CORONAL or SAGITALL orientation."
                       "Note that the reorientation is perform after the change of coordinate space if --space and --reorient"
                       "are given together.\n"
                       "INRIA / IRISA - VisAGeS Team",' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inputArg("i",
            "input",
            "input_filename",
            true,
            "",
            "Input image.",
            cmd);

    TCLAP::ValueArg<std::string> outputArg("o",
            "output",
            "output_filename",
            false,
            "",
            "Output image.",
            cmd);

    TCLAP::SwitchArg infoArg("I",
            "info",
            "Display info on the image or not",
            cmd);

    TCLAP::ValueArg<std::string> reorientArg("R",
            "reorient",
            "orientation",
            false,
            "",
            "Reorient the image in 'AXIAL' or 'CORONAL' or 'SAGITTAL' direction. [defalut: No reorientation]",
            cmd);
    TCLAP::ValueArg<std::string> spaceArg("s",
            "space",
            "space_filename",
            false,
            "",
            "Image used as space reference.",
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
    args.reorient = reorientArg.getValue(); args.space = spaceArg.getValue();

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
