#include <animaExpTensorImageFilter.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',"1.0");

    TCLAP::ValueArg<std::string> inArg("i","inputlist","Log-tensors image",true,"","input log-tensor image",cmd);
    TCLAP::ValueArg<std::string> resArg("o","outputfile","Result tensor image",true,"","result tensor image",cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("T","nb-threads","Number of threads (default: all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"Number of threads",cmd);
    TCLAP::SwitchArg scaleArg("S","scale","The log-tensors have their non-diagonal terms scaled",cmd,false);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    typedef anima::ExpTensorImageFilter <float> MainFilterType;

    typedef itk::ImageFileReader < MainFilterType::TInputImage > ImageReaderType;
    typedef itk::ImageFileWriter < MainFilterType::TOutputImage > ImageWriterType;

    MainFilterType::Pointer mainFilter = MainFilterType::New();

    mainFilter->SetNumberOfThreads(nbpArg.getValue());
    mainFilter->SetScaleNonDiagonal(scaleArg.isSet());

    ImageReaderType::Pointer tmpReader = ImageReaderType::New();
    tmpReader->SetFileName(inArg.getValue());
    tmpReader->Update();

    mainFilter->SetInput(tmpReader->GetOutput());

    mainFilter->Update();

    std::cout << "Writing result to : " << resArg.getValue() << std::endl;

    ImageWriterType::Pointer outWriter = ImageWriterType::New();
    outWriter->SetInput(mainFilter->GetOutput());
    outWriter->SetFileName(resArg.getValue());
    outWriter->SetUseCompression(true);

    outWriter->Update();

    return 0;
}
