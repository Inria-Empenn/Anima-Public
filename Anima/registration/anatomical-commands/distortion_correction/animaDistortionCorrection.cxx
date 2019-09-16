#include <animaDistortionCorrectionImageFilter.h>
#include <animaReadWriteFunctions.h>

#include <tclap/CmdLine.h>

#include <itkImage.h>

int main(int ac, const char** av)
{
    std::string descriptionMessage;
    descriptionMessage += "Compute a vector field in order to correct a distorted EPI\n";
    descriptionMessage += "INRIA / IRISA - VisAGeS/Empenn Team";

    TCLAP::CmdLine cmd(descriptionMessage, ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> forwardArg("f", "forwardImage", "Forward image", true, "", "Forward image", cmd);
    TCLAP::ValueArg<std::string> backwardArg("b", "backwardImage", "Backward image", true, "", "Backward image", cmd);
    TCLAP::ValueArg<unsigned int> distortionDirection("d", "dir", "Direction of distortion (0,1,2)", false, 1, "number of the direction of distortion", cmd);
    TCLAP::ValueArg<std::string> outArg("o", "outputVectorField", "Distortion vector field", true, "","Output vector field", cmd );

    TCLAP::ValueArg<double> sigmaArg("s", "sigma-smooth", "Sigma for gaussian smoothing of the transformation (in pixels, default: 2)",false,1,"gaussian smoothing sigma", cmd );
    TCLAP::ValueArg<unsigned int> nbpArg("T","numberofthreads","Number of threads to run on (default : all cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(ac,av);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    const unsigned int Dimension = 3;

    typedef double PixelType;

    typedef itk::Image <PixelType, Dimension> ImageType;
    typedef itk::VectorImage <PixelType, Dimension>  VectorFieldType;

    typedef anima::DistortionCorrectionImageFilter < ImageType > FilterType;

    ImageType::Pointer backwardImage = anima::readImage<ImageType>(backwardArg.getValue());
    backwardImage->DisconnectPipeline();

    double meanSpacing = 0;
    for (unsigned int i = 0;i < Dimension;++i)
        meanSpacing += backwardImage->GetSpacing()[i];

    meanSpacing /= Dimension;

    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(0, anima::readImage<ImageType>(forwardArg.getValue()));
    filter->SetInput(1, backwardImage);
    filter->SetDirection(distortionDirection.getValue());
    filter->SetFieldSmoothingSigma(sigmaArg.getValue() * meanSpacing);
    filter->SetNumberOfWorkUnits(nbpArg.getValue());
    filter->Update();

    anima::writeImage<VectorFieldType>(outArg.getValue(),filter->GetOutput());

    return EXIT_SUCCESS;
}
