#include <tclap/CmdLine.h>

#include <animaReadWriteFunctions.h>
#include <itkImage.h>
#include <itkVectorImage.h>

#include <animaSmoothingRecursiveYvvGaussianImageFilter.h>
#include <itkGradientRecursiveGaussianImageFilter.h>

template <class ImageType>
void performSmoothing (std::string &inStr, std::string &outStr, double sigma, unsigned int nThreads)
{
    typedef anima::SmoothingRecursiveYvvGaussianImageFilter <ImageType,ImageType> SmoothingFilterType;

    typename ImageType::Pointer currentImage = anima::readImage <ImageType> (inStr);
    currentImage->DisconnectPipeline();

    double meanSpacing = currentImage->GetSpacing()[0];
    for (unsigned int i = 1;i < ImageType::ImageDimension;++i)
        meanSpacing += currentImage->GetSpacing()[i];
    meanSpacing /= ImageType::ImageDimension;

    typename SmoothingFilterType::Pointer smoothFilter = SmoothingFilterType::New();

    smoothFilter->SetInput(currentImage);
    smoothFilter->SetNumberOfWorkUnits(nThreads);
    smoothFilter->SetSigma(sigma * meanSpacing);

    smoothFilter->Update();

    // Finally write the result
    anima::writeImage <ImageType> (outStr,smoothFilter->GetOutput());
}

template <class ImageType>
void performGradientEstimation (std::string &inStr, std::string &outStr, double sigma, unsigned int nThreads)
{
    typedef itk::Image <itk::Vector <typename ImageType::PixelType, ImageType::ImageDimension>, ImageType::ImageDimension> OutputImageType;
    typedef itk::GradientRecursiveGaussianImageFilter <ImageType,OutputImageType> GradientFilterType;

    typename ImageType::Pointer currentImage = anima::readImage <ImageType> (inStr);
    currentImage->DisconnectPipeline();

    double meanSpacing = currentImage->GetSpacing()[0];
    for (unsigned int i = 1;i < ImageType::ImageDimension;++i)
        meanSpacing += currentImage->GetSpacing()[i];
    meanSpacing /= ImageType::ImageDimension;

    typename GradientFilterType::Pointer gradientFilter = GradientFilterType::New();

    gradientFilter->SetInput(currentImage);
    gradientFilter->SetNumberOfWorkUnits(nThreads);
    gradientFilter->SetSigma(sigma * meanSpacing);

    gradientFilter->Update();

    // Finally write the result
    anima::writeImage <typename GradientFilterType::OutputImageType> (outStr,gradientFilter->GetOutput());
}

int main(int argc, char **argv)
{
    std::string descriptionMessage = "Performs Gaussian smoothing of an image\n";
    descriptionMessage += "INRIA / IRISA - VisAGeS/Empenn Team";

    TCLAP::CmdLine cmd(descriptionMessage, ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","input","Input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","output image",true,"","output image",cmd);

    TCLAP::ValueArg<double> sigmaArg("s","sigma","sigma value of Gaussian kernel",false,1.0,"sigma value",cmd);
    TCLAP::SwitchArg gradientArg("G","compute-gradient","Compute gradient instead of smoothed image",cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("T","numberofthreads","Number of threads to run on (default : all cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef itk::Image <double,3> ImageType;
    typedef itk::VectorImage <double,3> VectorImageType;

    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inArg.getValue().c_str(),
                                                                           itk::IOFileModeEnum::ReadMode);

    if(!imageIO)
    {
        std::cerr << "Itk could not find suitable IO factory for the input" << std::endl;
        return EXIT_FAILURE;
    }

    // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
    imageIO->SetFileName(inArg.getValue());
    imageIO->ReadImageInformation();

    if (!gradientArg.isSet())
    {
        bool vectorImage = (imageIO->GetNumberOfComponents() > 1);

        if (vectorImage)
            performSmoothing<VectorImageType>(inArg.getValue(),outArg.getValue(),sigmaArg.getValue(),nbpArg.getValue());
        else
            performSmoothing<ImageType>(inArg.getValue(),outArg.getValue(),sigmaArg.getValue(),nbpArg.getValue());
    }
    else
        performGradientEstimation <ImageType> (inArg.getValue(),outArg.getValue(),sigmaArg.getValue(),nbpArg.getValue());

    return EXIT_SUCCESS;
}
