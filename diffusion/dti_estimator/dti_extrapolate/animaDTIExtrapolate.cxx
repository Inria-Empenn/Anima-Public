#include <animaDTIExtrapolateImageFilter.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbe.h>

#include <tclap/CmdLine.h>

int main(int argc,  char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',"1.0");

    TCLAP::ValueArg<std::string> inArg("i","inputdti","DTI volume",true,"","DTI volume",cmd);
    TCLAP::ValueArg<std::string> resArg("o","output","Result DTI image",true,"","result DTI image",cmd);
    TCLAP::ValueArg<std::string> b0InArg("I","input-b0","Input B0 image",true,"","input B0 image",cmd);
    TCLAP::ValueArg<std::string> b0OutArg("O","output-b0","Resulting B0 image",true,"","result B0 image",cmd);

    TCLAP::ValueArg<unsigned int> nbFillsArg("n","numfills","Number of times filling is performed (default : 1)",false,1,"number of fills",cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default: all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    typedef anima::DTIExtrapolateImageFilter<float> FilterType;
    typedef FilterType::OutputImageType VectorImageType;
    typedef FilterType::OutputB0ImageType B0ImageType;

    typedef itk::ImageFileReader <VectorImageType> InputImageReaderType;
    typedef itk::ImageFileReader <B0ImageType> InputB0ImageReaderType;
    typedef itk::ImageFileWriter <B0ImageType> OutputB0ImageWriterType;
    typedef itk::ImageFileWriter <VectorImageType> VectorImageWriterType;

    InputImageReaderType::Pointer inputReader = InputImageReaderType::New();
    inputReader->SetFileName(inArg.getValue());

    inputReader->Update();
    VectorImageType::Pointer input = inputReader->GetOutput();

    InputB0ImageReaderType::Pointer inputB0Reader = InputB0ImageReaderType::New();
    inputB0Reader->SetFileName(b0InArg.getValue());

    inputB0Reader->Update();
    B0ImageType::Pointer inputB0 = inputB0Reader->GetOutput();

    itk::TimeProbe tmpTimer;

    tmpTimer.Start();

    for (unsigned int i = 0;i < nbFillsArg.getValue();++i)
    {
        FilterType::Pointer mainFilter = FilterType::New();

        mainFilter->SetInput(input);
        mainFilter->SetInitialEstimatedB0Image(inputB0);
        mainFilter->SetNumberOfThreads(nbpArg.getValue());

        try {
            mainFilter->Update();
        } catch (itk::ExceptionObject &e) {
            std::cerr << e << std::endl;
            exit(-1);
        }

        input = mainFilter->GetOutput();
        input->DisconnectPipeline();
        inputB0 = mainFilter->GetEstimatedB0Image();
    }

    tmpTimer.Stop();

    std::cout << "Extrapolation done in " << tmpTimer.GetTotal() << " s" << std::endl;

    std::cout << "Writing result to : " << resArg.getValue() << std::endl;

    VectorImageWriterType::Pointer vectorWriter = VectorImageWriterType::New();
    vectorWriter->SetInput(input);
    vectorWriter->SetFileName(resArg.getValue());
    vectorWriter->SetUseCompression(true);

    vectorWriter->Update();

    if (b0OutArg.getValue() != "")
    {
        OutputB0ImageWriterType::Pointer b0Writer = OutputB0ImageWriterType::New();
        b0Writer->SetInput(inputB0);
        b0Writer->SetFileName(b0OutArg.getValue());
        b0Writer->SetUseCompression(true);

        b0Writer->Update();
    }

    return 0;
}
