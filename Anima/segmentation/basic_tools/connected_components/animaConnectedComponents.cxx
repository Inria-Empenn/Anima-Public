#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>

#include <tclap/CmdLine.h>

#include <animaReadWriteFunctions.h>

struct arguments
{
    std::string input, output;
    unsigned int pthread;
    double min;
    bool full;
};

template <unsigned int Dimension>
void
connectedComponent(const arguments &args)
{
    typedef itk::Image <unsigned short,Dimension> ImageType;

    typedef itk::ConnectedComponentImageFilter <ImageType,ImageType> MainFilterType;
    typedef itk::RelabelComponentImageFilter <ImageType,ImageType> RelabelComponentFilterType;

    typename MainFilterType::Pointer mainFilter = MainFilterType::New();
    mainFilter->SetInput(anima::readImage<ImageType>(args.input));
    mainFilter->SetFullyConnected(args.full);
    if(args.pthread > 0)
        mainFilter->SetNumberOfThreads(args.pthread);
    mainFilter->Update();

    // Compute Image Spacing, 4th dimension is not physical but temporal
    typename ImageType::SpacingType spacing = mainFilter->GetInput(0)->GetSpacing();
    typename ImageType::SpacingValueType spacingTot = spacing[0];
    for (unsigned int i = 1;i < std::min(Dimension,(unsigned int)3);++i)
        spacingTot *= spacing[i];


    // Compute minsize in voxels
    double minSizeInVoxelD = args.min / spacingTot;
    double minSizeInVoxelD_floor = floor(minSizeInVoxelD);
    unsigned int minSizeInVoxel = static_cast<unsigned int>(minSizeInVoxelD_floor);
    minSizeInVoxel++; // to have strickly superior sizes
    double diff = minSizeInVoxelD-minSizeInVoxelD_floor;

    typename RelabelComponentFilterType::Pointer relabelFilter = RelabelComponentFilterType::New();
    relabelFilter->SetInput( mainFilter->GetOutput() );
    relabelFilter->SetMinimumObjectSize(minSizeInVoxel);
    if(args.pthread > 0)
        relabelFilter->SetNumberOfThreads(args.pthread);
    relabelFilter->Update();

    std::cout << "Original number of objects: " << relabelFilter->GetOriginalNumberOfObjects() <<std::endl;
    std::cout << "Total image spacing: " << spacingTot << std::endl;
    std::cout << "Connected components minimum size: " << args.min << " mm3 --> " << "process on " << minSizeInVoxel-1 << " voxel(s)" << std::endl;

    if (diff > 1.0e-6)
    {
        std::cout << "-- Warning: operation is not complete, " << (double)(minSizeInVoxel-1)*spacingTot << " mm3 is/are removed ("  << minSizeInVoxel-1  << " voxel(s)) "
                  << "instead of " << args.min << " mm3 (" << minSizeInVoxelD << " voxel(s))" << std::endl;
    }

    std::cout << "Number of objects after cleaning too small ones: " << relabelFilter->GetNumberOfObjects() << std::endl;
    std::cout << std::endl;

    anima::writeImage<ImageType>(args.output, relabelFilter->GetOutput());
}

void
retrieveDimension(const arguments &args, itk::ImageIOBase::Pointer imageIO)
{
    unsigned int nbDim = imageIO->GetNumberOfDimensions();

    switch(nbDim)
    {
        case 2:
            connectedComponent<2>(args);
            break;
        case 3:
            connectedComponent<3>(args);
            break;
        case 4:
            connectedComponent<4>(args);
            break;
        default:
            itk::ExceptionObject excp(__FILE__, __LINE__, "Number of dimension not supported.", ITK_LOCATION);
            throw excp;
    }
}

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","inputfile","Input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","Output image",true,"","output image",cmd);
    TCLAP::ValueArg<unsigned int> pArg("p","Thread","Number of thread",false,0,"Number of thread to use",cmd);

    TCLAP::ValueArg<double> minSizeArg("m","minsize","minimal component size in mm3",false,0,"minimal component size",cmd);
    TCLAP::SwitchArg fullConnectArg("F","full-connect","Use 26-connectivity instead of 6-connectivity",cmd,false);

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
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inArg.getValue().c_str(),
                                                                           itk::ImageIOFactory::ReadMode);
    if( !imageIO )
    {
        std::cerr << "Itk could not find suitable IO factory for the input" << std::endl;
        return EXIT_FAILURE;
    }

    // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
    imageIO->SetFileName(inArg.getValue());
    imageIO->ReadImageInformation();

    std::cout<<"\npreparing filter...\n";

    arguments args;

    args.input = inArg.getValue();
    args.output = outArg.getValue();
    args.min = minSizeArg.getValue();
    args.full= fullConnectArg.getValue();
    args.pthread = pArg.getValue();

    try
    {
        retrieveDimension(args, imageIO);
    }
    catch (itk::ExceptionObject & err)
    {
        std::cerr << "Itk cannot concatenate, be sure to use valid arguments..." << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }
}
