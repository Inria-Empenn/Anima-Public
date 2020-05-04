#include <tclap/CmdLine.h>

#include <itkImage.h>
#include <itkCommand.h>
#include <itkExtractImageFilter.h>

#include <animaReadWriteFunctions.h>

struct arguments
{
    std::vector<std::string> inputs;
    std::string base, output;
    double origin, spacing;
    unsigned int nbComponents;
};

template <class InputImageType, unsigned int InputDim, class OutputImageType >
void
concatenate(const arguments &args)
{
    unsigned int input = 0;

    typename OutputImageType::Pointer baseImg;
    if(args.base == "") // create args fake base from the first input.
    {
        std::cout << "No base given, direction of the created one will be guessed....\n" << std::endl;

        typename InputImageType::Pointer fimage = anima::readImage<InputImageType>(args.inputs[input]);
        ++input;

        typename OutputImageType::RegionType finalRegion;
        typename OutputImageType::SizeType finalSize;
        typename OutputImageType::IndexType finalIndex;
        typename OutputImageType::PointType finalOrigin;
        typename OutputImageType::SpacingType finalSpacing;
        typename OutputImageType::DirectionType finalDirection;

        for (unsigned int d = 0; d < InputDim; ++d)
        {

            finalIndex[d] = fimage->GetLargestPossibleRegion().GetIndex()[d];
            finalSize[d] = fimage->GetLargestPossibleRegion().GetSize()[d];
            finalOrigin[d] = fimage->GetOrigin()[d];
            finalSpacing[d] = fimage->GetSpacing()[d];
            for(unsigned int i = 0; i < InputDim; ++i)
                finalDirection[d][i] = fimage->GetDirection()[d][i];
        }
        finalIndex[InputDim] = 0;
        finalSize[InputDim] = 1;
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

        baseImg = OutputImageType::New();
        baseImg->Initialize();
        baseImg->SetRegions(finalRegion);
        baseImg->SetOrigin(finalOrigin);
        baseImg->SetSpacing(finalSpacing);
        baseImg->SetDirection(finalDirection);
        baseImg->SetNumberOfComponentsPerPixel(args.nbComponents);
        baseImg->Allocate();

        typedef itk::ImageRegionIterator <OutputImageType> FillIteratorType;
        FillIteratorType fillItr(baseImg, baseImg->GetLargestPossibleRegion());
        typedef itk::ImageRegionConstIterator<InputImageType> SourceIteratorType;
        SourceIteratorType srcItr(fimage, fimage->GetLargestPossibleRegion());
        while(!srcItr.IsAtEnd())
        {
            fillItr.Set(srcItr.Get());
            ++fillItr; ++srcItr;
        }

        std::cout << "input " << input << " / " << args.inputs.size() << " concatenated..." << std::endl;
    }
    else // use given basis.
        baseImg = anima::readImage<OutputImageType>(args.base);

    // allocate output
    typename OutputImageType::Pointer outputImg;
    typename OutputImageType::RegionType finalRegion;
    typename OutputImageType::SizeType finalSize;
    typename OutputImageType::IndexType finalIndex;
    typename OutputImageType::PointType finalOrigin;
    typename OutputImageType::SpacingType finalSpacing;
    typename OutputImageType::DirectionType finalDirection;

    for (unsigned int d = 0; d < InputDim; ++d)
    {
        finalIndex[d] = baseImg->GetLargestPossibleRegion().GetIndex()[d];
        finalSize[d] = baseImg->GetLargestPossibleRegion().GetSize()[d];
        finalOrigin[d] = baseImg->GetOrigin()[d];
        finalSpacing[d] = baseImg->GetSpacing()[d];

        for(unsigned int i = 0;i < InputDim;++i)
            finalDirection[d][i] = baseImg->GetDirection()[d][i];
    }

    for(unsigned int i = 0;i < InputDim + 1;++i)
    {
        finalDirection[InputDim][i] = baseImg->GetDirection()[InputDim][i];
        finalDirection[i][InputDim] = baseImg->GetDirection()[i][InputDim];
    }

    finalIndex[InputDim] = baseImg->GetLargestPossibleRegion().GetIndex()[InputDim];
    finalSize[InputDim] = baseImg->GetLargestPossibleRegion().GetSize()[InputDim] + args.inputs.size() - input;
    finalOrigin[InputDim] = baseImg->GetOrigin()[InputDim];
    finalSpacing[InputDim] = baseImg->GetSpacing()[InputDim];

    finalRegion.SetIndex(finalIndex);
    finalRegion.SetSize(finalSize);

    outputImg = OutputImageType::New();
    outputImg->Initialize();
    outputImg->SetRegions(finalRegion);
    outputImg->SetOrigin(finalOrigin);
    outputImg->SetSpacing(finalSpacing);
    outputImg->SetDirection(finalDirection);
    outputImg->SetNumberOfComponentsPerPixel(args.nbComponents);
    outputImg->Allocate();

    // fill with the base:
    typedef itk::ImageRegionIterator <OutputImageType> FillIteratorType;
    FillIteratorType fillItr(outputImg, outputImg->GetLargestPossibleRegion());
    typedef itk::ImageRegionConstIterator<OutputImageType> BaseIteratorType;
    BaseIteratorType baseItr(baseImg, baseImg->GetLargestPossibleRegion());
    while(!baseItr.IsAtEnd())
    {
        fillItr.Set(baseItr.Get());
        ++fillItr; ++baseItr;
    }

    //fill with the inputs:
    for(unsigned int i = input; i < args.inputs.size();++i)
    {
        typename InputImageType::Pointer inputImg = anima::readImage<InputImageType>(args.inputs[i]);

        typedef itk::ImageRegionConstIterator<InputImageType> InputIteratorType;
        InputIteratorType inputItr(inputImg, inputImg->GetLargestPossibleRegion());
        while(!inputItr.IsAtEnd())
        {
            fillItr.Set(inputItr.Get());
            ++fillItr; ++inputItr;
        }

        std::cout << "input " << i+1 << " / " << args.inputs.size() << " concatenated..." << std::endl;
    }

    // write res
    anima::writeImage<OutputImageType>(args.output, outputImg);
}

template <class ComponentType, unsigned int InputDim>
void
retrieveNbComponent(arguments &args, itk::ImageIOBase::Pointer imageIO)
{
    args.nbComponents = imageIO->GetNumberOfComponents();
    switch(args.nbComponents)
    {
    case 1:
        typedef itk::Image<ComponentType, InputDim> InputImageType;
        typedef itk::Image<ComponentType, InputDim + 1> OutputImageType;
        concatenate <InputImageType, InputDim, OutputImageType >(args);
        break;
    case 3:
    case 6:
        typedef itk::VectorImage<ComponentType, InputDim> InputVectorImageType;
        typedef itk::VectorImage<ComponentType, InputDim + 1> OutputVectorImageType;
        concatenate < InputVectorImageType, InputDim, OutputVectorImageType >(args);
    break;
    default:
        itk::ExceptionObject excp(__FILE__, __LINE__, "Number of component not supported.", ITK_LOCATION);
        throw excp;
    }
}

template <class ComponentType >
void
retrieveNbDimensions(arguments &args, itk::ImageIOBase::Pointer imageIO)
{

    unsigned int nbDim = imageIO->GetNumberOfDimensions();

    switch(nbDim)
    {
    case 2:
        retrieveNbComponent<ComponentType, 2>(args, imageIO);
        break;
    case 3:
        retrieveNbComponent<ComponentType, 3>(args, imageIO);
        break;
    default:
        itk::ExceptionObject excp(__FILE__, __LINE__, "Number of dimension not supported.", ITK_LOCATION);
        throw excp;
    }
}

void
retrieveComponentType(arguments &args, itk::ImageIOBase::Pointer imageIO)
{
    switch (imageIO->GetComponentType())
    {
    case itk::ImageIOBase::UCHAR:
    case itk::ImageIOBase::USHORT:
    case itk::ImageIOBase::UINT:
    case itk::ImageIOBase::ULONG:
        retrieveNbDimensions<unsigned int>(args, imageIO);
        break;
    case itk::ImageIOBase::CHAR:
    case itk::ImageIOBase::SHORT:
    case itk::ImageIOBase::INT:
    case itk::ImageIOBase::LONG:
        retrieveNbDimensions<int>(args, imageIO);
        break;
    case itk::ImageIOBase::FLOAT:
        retrieveNbDimensions<double>(args, imageIO);
        break;
    case itk::ImageIOBase::DOUBLE:
        retrieveNbDimensions<double>(args, imageIO);
        break;
    default:
        itk::ExceptionObject excp(__FILE__, __LINE__, "Component type not supported.", ITK_LOCATION);
        throw excp;
    }
}


int main(int ac, const char** av)
{

    TCLAP::CmdLine cmd("Concatenate args serie of image between them.\nYou can give several -i input to concatenate.\n\
A -b base image is facultative, if it is given, all the inputs will be add at the end of the last dim of this base.\n\
if args base is given the geometry of the base is used otherwise it use the geometry of the first given input.\n\
with the origin -O(def = 0) and spacing -s(def = 1) passed as arguments.\
It's your responsability to give args set of well formed input.\n\n\
Example: the arguments -b 4x4x4x1 and -i 4x4x4 -i 4x4x4 will result on an outpu 4x4x4x3.\n\
INRIA / IRISA - VisAGeS/Empenn Team",
            ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> baseArg("b",
            "base",
            "Inputs will be concatenate with this base image",
            false,
            "",
            "Base image.",
            cmd);

    TCLAP::MultiArg<std::string> inputArg("i",
            "input",
            "Input images to concatenate",
            true,
            "Inputs",
            cmd);

    TCLAP::ValueArg<std::string> outputArg("o",
            "output",
            "Output image",
            true,
            "",
            "Output",
            cmd);
    TCLAP::ValueArg<double> originArg("O",
            "origin",
            "Origin of the last dimension of the concatenated image. ignore if args base has been given",
            false,
            0,
            "Origin",
            cmd);
    TCLAP::ValueArg<double> spacingArg("s",
            "spacing",
            "Spacing of the last dimension of the concatenated image. ignore if args base has been given",
            false,
            1,
            "Spacing",
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

    arguments args;
    args.inputs = inputArg.getValue();

    // Find out the type of the image in file
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(args.inputs[0].c_str(),
                                                                           itk::ImageIOFactory::ReadMode);

    if( !imageIO )
    {
        std::ifstream fileIn(args.inputs[0].c_str());
        std::vector <std::string> listInputNames;

        if (!fileIn.is_open())
        {
            std::cerr << "Unable to read file: " << args.inputs[0] << std::endl;
            return EXIT_FAILURE;
        }

        while (!fileIn.eof())
        {
            char tmpStr[2048];
            fileIn.getline(tmpStr,2048);

            if (strcmp(tmpStr,"") == 0)
                continue;

            listInputNames.push_back(tmpStr);
        }

        fileIn.close();

        if (listInputNames.size() > 0)
            imageIO = itk::ImageIOFactory::CreateImageIO(listInputNames[0].c_str(),itk::ImageIOFactory::ReadMode);

        if( !imageIO )
        {
            std::cerr << "Itk could not find suitable IO factory for the input" << std::endl;
            return EXIT_FAILURE;
        }

        args.inputs = listInputNames;
    }

    // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
    imageIO->SetFileName(args.inputs[0]);
    imageIO->ReadImageInformation();

    std::cout<<"\npreparing filter...\n";

    args.output = outputArg.getValue();
    args.base = baseArg.getValue();
    args.origin = originArg.getValue();
    args.spacing = spacingArg.getValue();

    try
    {
        retrieveComponentType(args, imageIO);
    }
    catch ( itk::ExceptionObject & err )
    {
        std::cerr << "Itk cannot concatenate, be sure to use valid arguments..." << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
