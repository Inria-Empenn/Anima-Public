#include <tclap/CmdLine.h>

#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkCommand.h>
#include <itkExtractImageFilter.h>

#include <animaReadWriteFunctions.h>
#include <animaRetrieveImageTypeMacros.h>

struct arguments
{
    std::string input, output;
    double origin, spacing;
};

template <class ComponentType, unsigned int InputDim>
void
collapse(itk::ImageIOBase::Pointer imageIO, const arguments &args)
{

    typedef itk::VectorImage<ComponentType, InputDim> InputImageType;
    typedef itk::Image<ComponentType, InputDim+1> OutputImageType;

    typename InputImageType::Pointer inputImg = anima::readImage<InputImageType>(args.input);
    unsigned int nbComp = imageIO->GetNumberOfComponents();

    std::cout << "number of dimension : " << InputDim << std::endl;
    std::cout << "number of components : " << nbComp << std::endl;


    typename OutputImageType::RegionType finalRegion;
    typename OutputImageType::SizeType finalSize;
    typename OutputImageType::IndexType finalIndex;
    typename OutputImageType::PointType finalOrigin;
    typename OutputImageType::SpacingType finalSpacing;
    typename OutputImageType::DirectionType finalDirection;

    for (unsigned int d = 0; d < InputDim; ++d)
    {

        finalIndex[d] = inputImg->GetLargestPossibleRegion().GetIndex()[d];
        finalSize[d] = inputImg->GetLargestPossibleRegion().GetSize()[d];
        finalOrigin[d] = inputImg->GetOrigin()[d];
        finalSpacing[d] = inputImg->GetSpacing()[d];
        for(unsigned int i = 0; i < InputDim; ++i)
            finalDirection[d][i] = inputImg->GetDirection()[d][i];
    }
    finalIndex[InputDim] = 0;
    finalSize[InputDim] = nbComp;
    finalOrigin[InputDim] = args.origin;
    finalSpacing[InputDim] = args.spacing;

    for(unsigned int i = 0; i < InputDim; ++i)
    {
        finalDirection[InputDim][i] = 0;
        finalDirection[i][InputDim] = 0;
    }
    finalDirection[InputDim][InputDim] = 1;

    finalRegion.SetIndex(finalIndex);
    finalRegion.SetSize(finalSize);

    typename OutputImageType::Pointer outputImg = OutputImageType::New();
    outputImg->Initialize();
    outputImg->SetRegions(finalRegion);
    outputImg->SetOrigin(finalOrigin);
    outputImg->SetSpacing(finalSpacing);
    outputImg->SetDirection(finalDirection);
    outputImg->Allocate();

    typedef itk::ImageRegionIterator <OutputImageType> FillIteratorType;
    FillIteratorType fillItr(outputImg, outputImg->GetLargestPossibleRegion());
    typedef itk::ImageRegionConstIterator<InputImageType> SourceIteratorType;
    SourceIteratorType srcItr(inputImg, inputImg->GetLargestPossibleRegion());
    int idx = 0;
    while(!fillItr.IsAtEnd())
    {
        if(srcItr.IsAtEnd())
        {
            idx++;
            srcItr.GoToBegin();
        }

        fillItr.Set(srcItr.Get()[idx]);
        ++fillItr; ++srcItr;
    }

    anima::writeImage<OutputImageType>(args.output, outputImg);
}

template <class ComponentType >
void
retrieveNbDimensions(itk::ImageIOBase::Pointer imageIO, const arguments &args)
{
    ANIMA_RETRIEVE_NUMBER_OF_DIMENSIONS(imageIO, ComponentType, collapse, imageIO, args)
}


int main(int ac, const char** av)
{

    TCLAP::CmdLine cmd("Collapse a vector image into a scalar image of dim n+1.\n\
The last dimension of the output has same size as the component of the input image.\n\
You can give the spacing and the origin of the added dimension, default values are 0 for the origin and 1 for the spacing\n\
Example: input image is 4x4x4 with a component size of 32, the output will be 4x4x4x32.\n\
INRIA / IRISA - VisAGeS/Empenn Team",
            ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inputArg("i",
            "input",
            "Input images to collapse",
            true,
            "",
            "Inputs",
            cmd);

    TCLAP::ValueArg<std::string> outputArg("o",
            "output",
            "Output image",
            true,
            "",
            "Output",
            cmd);
    TCLAP::ValueArg<double> origArg("O",
            "origin",
            "Origin of the added Dimension image",
            false,
            0,
            "Origin",
            cmd);
    TCLAP::ValueArg<double> spacingArg("s",
            "spacing",
            "Spacing of the added Dimension image",
            false,
            1,
            "Spacing ",
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

    if(imageIO->GetNumberOfComponents() < 2)
    {
        std::cerr << "Input have to have component of size > 1. Nothing to do.";
        return EXIT_SUCCESS;
    }

    std::cout<<"\npreparing filter...\n";

    arguments args;

    args.input = inputArg.getValue();
    args.output = outputArg.getValue();
    args.origin = origArg.getValue();
    args.spacing = spacingArg.getValue();


    try
    {
        ANIMA_RETRIEVE_COMPONENT_TYPE(imageIO, retrieveNbDimensions, imageIO, args);
    }
    catch ( itk::ExceptionObject & err )
    {
        std::cerr << "Itk cannot collapse, be sure to use valid input..." << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
