#include <tclap/CmdLine.h>

#include <itkHistogramMatchingImageFilter.h>

#include <animaReadWriteFunctions.h>
#include <animaRetrieveImageTypeMacros.h>

struct arguments
{
    std::string reference, moving, output;
    unsigned int nlevels, npoints, nthreads;
};

template <class ComponentType, unsigned int InputDim>
void
standardize(itk::ImageIOBase::Pointer imageIO, const arguments &args)
{
    
    typedef itk::Image<ComponentType, InputDim> ImageType;
    typedef itk::HistogramMatchingImageFilter<ImageType,ImageType> HistogramMatchingFilterType;
    
    typename HistogramMatchingFilterType::Pointer mainFilter = HistogramMatchingFilterType::New();
    mainFilter->SetReferenceImage(anima::readImage<ImageType>(args.reference));
    mainFilter->SetInput(anima::readImage<ImageType>(args.moving));
    mainFilter->SetNumberOfHistogramLevels(args.nlevels);
    mainFilter->SetNumberOfMatchPoints(args.npoints);
    mainFilter->SetNumberOfThreads(args.nthreads);
    mainFilter->Update();
    
    anima::writeImage<ImageType>(args.output, mainFilter->GetOutput());
}

template <class ComponentType>
void
retrieveNbDimensions(itk::ImageIOBase::Pointer imageIO, const arguments &args)
{
    ANIMA_RETRIEVE_NUMBER_OF_DIMENSIONS(imageIO, ComponentType, standardize, imageIO, args)
}

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> refArg("r","referencefile","Reference image",true,"","reference image",cmd);
    TCLAP::ValueArg<std::string> movArg("m","movingfile","Moving image",true,"","moving image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","output image",true,"","output image",cmd);
    
    TCLAP::ValueArg<unsigned int> nlvlArg("","nhl","Number of histogram levels (default: 100)",false,100,"number of histogram levels",cmd);
    TCLAP::ValueArg<unsigned int> nptsArg("","nmp","Number of match points (default: 15)",false,15,"number of match points",cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default: all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    // Find out the type of the image in file
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(refArg.getValue().c_str(),itk::ImageIOFactory::ReadMode);

    if (!imageIO)
    {
        std::cerr << "Itk could not find a suitable IO factory for the input" << std::endl;
        return EXIT_FAILURE;
    }

    // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
    imageIO->SetFileName(refArg.getValue());
    imageIO->ReadImageInformation();
    
    arguments args;
    
    args.reference = refArg.getValue();
    args.moving = movArg.getValue();
    args.output = outArg.getValue();
    args.nlevels = nlvlArg.getValue();
    args.npoints = nptsArg.getValue();
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
}
