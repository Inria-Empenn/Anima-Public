#include <itkSignedMaurerDistanceMapImageFilter.h>

#include <animaReadWriteFunctions.h>
#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inArg("i","inputfile","Input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","Output image",true,"","output image",cmd);

    TCLAP::SwitchArg invArg("I","positive-outside","Positive distances will be outside the object (default: not activated)",cmd,false);
    
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef itk::Image <unsigned short,3> ImageType;
    typedef itk::Image <float,3> FloatImageType;

    typedef itk::SignedMaurerDistanceMapImageFilter <ImageType,FloatImageType> MainFilterType;

    MainFilterType::Pointer mainFilter = MainFilterType::New();
    mainFilter->SetInput(anima::readImage <ImageType> (inArg.getValue()));
    mainFilter->SetSquaredDistance(false);
    mainFilter->SetBackgroundValue(0);
    
    if (!invArg.isSet())
        mainFilter->InsideIsPositiveOn();
    else
        mainFilter->InsideIsPositiveOff();

    mainFilter->UseImageSpacingOn();
    mainFilter->Update();
    
    anima::writeImage <FloatImageType> (outArg.getValue(),mainFilter->GetOutput());
    return EXIT_SUCCESS;
}
