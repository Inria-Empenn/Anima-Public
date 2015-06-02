#include <animaFDRCorrectImageFilter.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',"1.0");

    TCLAP::ValueArg<std::string> inArg("i","input","Non corrected P-value image",true,"","Non corrected P-value image",cmd);
    TCLAP::ValueArg<std::string> resArg("o","output","FDR corrected output image",true,"","fdr corrected output image",cmd);
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

    typedef itk::ImageFileReader < MainFilterType::TInputImage > ImageReaderType;
    typedef itk::ImageFileReader < MainFilterType::MaskImageType > MaskImageReaderType;
    typedef itk::ImageFileWriter < MainFilterType::TOutputImage > ImageWriterType;

    MainFilterType::Pointer mainFilter = MainFilterType::New();

    ImageReaderType::Pointer tmpReader = ImageReaderType::New();
    tmpReader->SetFileName(inArg.getValue());
    tmpReader->Update();
    mainFilter->SetInput(tmpReader->GetOutput());

    if (maskArg.getValue() != "")
    {
        MaskImageReaderType::Pointer maskReader = MaskImageReaderType::New();
        maskReader->SetFileName(maskArg.getValue());
        maskReader->Update();
        mainFilter->SetMaskImage(maskReader->GetOutput());
    }

    mainFilter->Update();

    std::cout << "Writing result to : " << resArg.getValue() << std::endl;

    ImageWriterType::Pointer outWriter = ImageWriterType::New();
    outWriter->SetInput(mainFilter->GetOutput());
    outWriter->SetFileName(resArg.getValue());
    outWriter->SetUseCompression(true);

    outWriter->Update();

    return 0;
}
