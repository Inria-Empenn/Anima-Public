#include <tclap/CmdLine.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkCommand.h>
#include <itkExtractImageFilter.h>

#include <animaReadWriteFunctions.h>
#include <animaRetrieveImageTypeMacros.h>


//Update progression of the process
void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
    std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
}

struct arguments
{
    int xindex, xsize, yindex, ysize, zindex, zsize, tindex, tsize;
    std::string input, output;
};


template <class InputImageType, unsigned int OutputDimension>
void
extract(const arguments &args)
{
    typename InputImageType::Pointer input = anima::readImage<InputImageType>(args.input);

    typename InputImageType::SizeType inputSize;
    inputSize = input->GetLargestPossibleRegion().GetSize();

    typename InputImageType::RegionType extractRegion;
    typename InputImageType::SizeType extractSize;
    typename InputImageType::IndexType extractIndex;

    int indexes[4] = {args.xindex, args.yindex, args.zindex, args.tindex};
    int sizes[4] = {args.xsize, args.ysize, args.zsize, args.tsize};

    for (unsigned int d = 0; d < InputImageType::ImageDimension; ++d)
    {
        extractIndex[d] = indexes[d];
        if(sizes[d] == -1 || indexes[d] + sizes[d] > inputSize[d])
            extractSize[d] = inputSize[d] - extractIndex[d];
        else
            extractSize[d] = sizes[d];
    }

    extractRegion.SetIndex(extractIndex);
    extractRegion.SetSize(extractSize);

    std::cout<< "Input will be crop using ROI of dimensions:" << extractRegion << std::endl;

    if(input->GetNumberOfComponentsPerPixel() == 1)
    {
        typedef itk::Image<typename InputImageType::PixelType, OutputDimension> OutputImageType;

        typedef itk::ExtractImageFilter<InputImageType, OutputImageType> ExtractFilterType;
        typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
        extractFilter->SetExtractionRegion(extractRegion);
        extractFilter->SetDirectionCollapseToGuess();
        extractFilter->SetInput(input);

        itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
        callback->SetCallback(eventCallback);
        extractFilter->AddObserver(itk::ProgressEvent(), callback);
        extractFilter->Update();

        typename OutputImageType::Pointer output = extractFilter->GetOutput();

        std::cout << "\n\nOutput dimensions: " << output->GetLargestPossibleRegion();
        anima::writeImage<OutputImageType>(args.output, output);
    }
    else
    {
        typedef itk::VectorImage<typename InputImageType::InternalPixelType, OutputDimension> OutputImageType;

        typedef itk::ExtractImageFilter<InputImageType, OutputImageType> ExtractFilterType;
        typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
        extractFilter->SetExtractionRegion(extractRegion);
        extractFilter->SetDirectionCollapseToGuess();
        extractFilter->SetInput(input);

        itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
        callback->SetCallback(eventCallback);
        extractFilter->AddObserver(itk::ProgressEvent(), callback);
        extractFilter->Update();

        typename OutputImageType::Pointer output = extractFilter->GetOutput();

        std::cout << "\n\nOutput dimensions: " << output->GetLargestPossibleRegion();
        anima::writeImage<OutputImageType>(args.output, output);
    }


}

template <class InputImageType>
void
evaluateOutputType(const arguments &args)
{
    unsigned int imageDim = InputImageType::ImageDimension;
    switch(imageDim)
    {
    case 2:
        {
            unsigned int outputDim = 2;
            outputDim -= (args.xsize == 0)? 1 :0;
            outputDim -= (args.ysize == 0)? 1 :0;
            switch(outputDim)
            {
            case 1:
                extract<InputImageType, 1>(args);
                break;
            case 2:
                extract<InputImageType, 2>(args);
                break;
            default:
                std::string msg = "Number of collapsed dimension not supported.";
                itk::ExceptionObject excp(__FILE__, __LINE__,msg , ITK_LOCATION);
                throw excp;
            }
            break;
        }
        case 3:
        {
            unsigned int outputDim = 3;
            outputDim -= (args.xsize == 0)? 1 :0;
            outputDim -= (args.ysize == 0)? 1 :0;
            outputDim -= (args.zsize == 0)? 1 :0;
            switch(outputDim)
            {
            case 1:
                extract<InputImageType, 1>(args);
                break;
            case 2:
                extract<InputImageType, 2>(args);
                break;
            case 3:
                extract<InputImageType, 3>(args);
                break;
            default:
                std::string msg = "Number of collapsed dimension not supported.";
                itk::ExceptionObject excp(__FILE__, __LINE__,msg , ITK_LOCATION);
                throw excp;
            }
            break;
        }
        case 4:
        {
            unsigned int outputDim = 4;
            outputDim -= (args.xsize == 0)? 1 :0;
            outputDim -= (args.ysize == 0)? 1 :0;
            outputDim -= (args.zsize == 0)? 1 :0;
            outputDim -= (args.tsize == 0)? 1 :0;
            switch(outputDim)
            {
            case 1:
                extract<InputImageType, 1>(args);
                break;
            case 2:
                extract<InputImageType, 2>(args);
                break;
            case 3:
                extract<InputImageType, 3>(args);
                break;
            case 4:
                extract<InputImageType, 4>(args);
                break;
            default:
                std::string msg = "Number of collapsed dimension not supported.";
                itk::ExceptionObject excp(__FILE__, __LINE__,msg , ITK_LOCATION);
                throw excp;
            }
            break;
        }
        default:
            std::string msg = "Number of collapsed dimension not supported.";
            itk::ExceptionObject excp(__FILE__, __LINE__,msg , ITK_LOCATION);
            throw excp;
    }
}

template <class ComponentType, int dimension>
void
checkIfComponentsAreVectors(itk::ImageIOBase::Pointer imageIO, const arguments &args)
{
    ANIMA_CHECK_IF_COMPONENTS_ARE_VECTORS(imageIO, ComponentType, dimension, evaluateOutputType, args)
}


template <class ComponentType >
void
retrieveNbDimensions(itk::ImageIOBase::Pointer imageIO, const arguments &args)
{
    ANIMA_RETRIEVE_NUMBER_OF_DIMENSIONS(imageIO, ComponentType, checkIfComponentsAreVectors, imageIO, args)
}

int main(int ac, const char** av)
{

    TCLAP::CmdLine cmd("The animaCropImage uses an itkExtractImage filter to crop "
                       "an image given as input.\n"
                       "The lower case arguments(x<xindex>, y<yindex>, z<zindex>, t<tindex>)"
                       " are the starting indexes of the input region to keep. The default value is 0\n"
                       "The upper case arguments(X<xsize>, Y<ysize>, Z<zsize>, T<tsize>) are the sizes of "
                       "the input region to keep. The default value is the largest possible sizes given the "
                       "corresponding indexes.\nIf you give args size of zero the corresponding dimension will "
                       "be collapsed.\nExample: for args a 4D image 4x4x4x4 the arguments :\n --xindex 1"
                       " --zindex 1 --zsize 2 --tindex 3 --tsize 0\n will result on an image 3x4x2\n"
                       "Where the x dim corresponds to [1,2,3] of the input, y[0,3], zindex[1,2] and tindex is "
                       "collapsed, only the last sequence has been kept.",
                       ' ',
                       "0.0");

    TCLAP::ValueArg<std::string> inputArg("i",
            "input",
            "Input image to crop",
            true,
            "",
            "Input image to crop",
            cmd);

    TCLAP::ValueArg<std::string> outputArg("o",
            "output",
            "Output cropped image",
            true,
            "",
            "Output cropped image",
            cmd);

    TCLAP::ValueArg<unsigned int> xArg("x",
            "xindex",
            "The resulting croped image will go from xindex to xsize along the xindex axis.",
            false,
            0,
            "Start of ROI for the xindex dimension",
            cmd);

    TCLAP::ValueArg<unsigned int> XArg("X",
            "xsize",
            "The resulting croped image will go from xindex to xsize along the xindex axis. If 0 the dimension is collapsed.",
            false,
            -1,
            "Size of ROI for the xindex dimension",
            cmd);

    TCLAP::ValueArg<unsigned int> yArg("y",
            "yindex",
            "The resulting croped image will go from yindex to ysize along the yindex axis.",
            false,
            0,
            "Start of ROI for the yindex dimension",
            cmd);

    TCLAP::ValueArg<unsigned int> YArg("Y",
            "ysize",
            "The resulting croped image will go from yindex to ysize along the yindex axis. If 0 the dimension is collapsed.",
            false,
            -1,
            "Size of ROI for the yindex dimension",
            cmd);
    TCLAP::ValueArg<unsigned int> zArg("z",
            "zindex",
            "The resulting croped image will go from yindex to ysize along the zindex axis.",
            false,
            0,
            "Start of ROI for the zindex dimension",
            cmd);

    TCLAP::ValueArg<unsigned int> ZArg("Z",
            "zsize",
            "The resulting croped image will go from zindex to zsize along the zindex axis. If 0 the dimension is collapsed.",
            false,
            -1,
            "Size of ROI for the zindex dimension",
            cmd);
    TCLAP::ValueArg<unsigned int> tArg("t",
            "tindex",
            "The resulting croped image will go from tindex to tsize along the tindex axis.",
            false,
            0,
            "Start of ROI for the tindex dimension",
            cmd);

    TCLAP::ValueArg<unsigned int> TArg("T",
            "tsize",
            "The resulting croped image will go from tindex to tsize along the tindex axis. If 0 the dimension is collapsed.",
            false,
            -1,
            "Size of ROI for the tindex dimension",
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

    std::cout<<"\npreparing filter...\n";

    arguments args;
    args.xindex = xArg.getValue(); args.xsize = XArg.getValue();
    args.yindex = yArg.getValue(); args.ysize = YArg.getValue();
    args.zindex = zArg.getValue(); args.zsize = ZArg.getValue();
    args.tindex = tArg.getValue(); args.tsize = TArg.getValue();
    args.input = inputArg.getValue(); args.output = outputArg.getValue();

    try
    {
        ANIMA_RETRIEVE_COMPONENT_TYPE(imageIO, retrieveNbDimensions, imageIO, args);
    }
    catch ( itk::ExceptionObject & err )
    {
        std::cerr << "Itk cannot extract, be sure to use valid arguments..." << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
