#include <tclap/CmdLine.h>

#include <itkHistogramMatchingImageFilter.h>
#include <animaReadWriteFunctions.h>

using namespace std;

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> refArg("r","referencefile","Reference image",true,"","reference image",cmd);
    TCLAP::ValueArg<std::string> movArg("m","movingfile","Moving image",true,"","moving image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","output image",true,"","output image",cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default : all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

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

    if( !imageIO )
    {
        std::cerr << "Itk could not find a suitable IO factory for the input" << std::endl;
        return EXIT_FAILURE;
    }

    // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
    imageIO->SetFileName(refArg.getValue());
    imageIO->ReadImageInformation();

    const unsigned int nbDimension = imageIO->GetNumberOfDimensions();

    switch (nbDimension)
    {
        case 4:
        {
            typedef itk::Image<float,4> ImageType;
            typedef itk::HistogramMatchingImageFilter<ImageType,ImageType> HistogramMatchingFilterType;
            
            HistogramMatchingFilterType::Pointer mainFilter = HistogramMatchingFilterType::New();
            mainFilter->SetReferenceImage(anima::readImage<ImageType>(refArg.getValue()));
            mainFilter->SetInput(anima::readImage<ImageType>(movArg.getValue()));
            mainFilter->SetNumberOfHistogramLevels(100);
            mainFilter->SetNumberOfMatchPoints(15);
            mainFilter->SetNumberOfThreads(nbpArg.getValue());
            mainFilter->Update();

            anima::writeImage<ImageType>(outArg.getValue(),mainFilter->GetOutput());

            break;
        }

        case 3:
        default:
        {
            typedef itk::Image<float,3> ImageType;
            typedef itk::HistogramMatchingImageFilter<ImageType,ImageType> HistogramMatchingFilterType;
            
            HistogramMatchingFilterType::Pointer mainFilter = HistogramMatchingFilterType::New();
            mainFilter->SetReferenceImage(anima::readImage<ImageType>(refArg.getValue()));
            mainFilter->SetInput(anima::readImage<ImageType>(movArg.getValue()));
            mainFilter->SetNumberOfHistogramLevels(100);
            mainFilter->SetNumberOfMatchPoints(15);
            mainFilter->SetNumberOfThreads(nbpArg.getValue());
            mainFilter->Update();
            
            anima::writeImage<ImageType>(outArg.getValue(),mainFilter->GetOutput());
            
            break;
        }
    }

    return EXIT_SUCCESS;
}
