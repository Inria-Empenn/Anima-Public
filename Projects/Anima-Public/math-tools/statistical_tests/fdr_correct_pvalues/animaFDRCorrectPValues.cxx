#include <animaFDRCorrectImageFilter.h>
#include <animaReadWriteFunctions.h>

#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',"1.0");

    TCLAP::ValueArg<std::string> inArg("i","input","Non corrected P-value image",true,"","Non corrected P-value image",cmd);
    TCLAP::ValueArg<std::string> resArg("o","output","FDR thresholded output image at q",true,"","FDR corrected output image at q",cmd);
    TCLAP::ValueArg<double> qArg("q","q-val","FDR q value",true,0.05,"FDR q value",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","mask","Mask image (default: all pixels are in mask)",false,"","Mask image",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    typedef anima::FDRCorrectImageFilter<double> MainFilterType;

    MainFilterType::Pointer mainFilter = MainFilterType::New();

    mainFilter->SetInput(anima::readImage<MainFilterType::TInputImage> (inArg.getValue()));
    mainFilter->SetQValue(qArg.getValue());

    if (maskArg.getValue() != "")
        mainFilter->SetMaskImage(anima::readImage <MainFilterType::MaskImageType> (maskArg.getValue()));

    mainFilter->Update();

    std::cout << "Writing result to : " << resArg.getValue() << std::endl;

    anima::writeImage<MainFilterType::TOutputImage>(resArg.getValue(), mainFilter->GetOutput());

    return 0;
}
