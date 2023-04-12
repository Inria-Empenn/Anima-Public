#include <animaSampleImageFromDistributionImageFilter.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
	
    TCLAP::ValueArg<std::string> meanArg("m","mean","Average of the distribution",true,"","Average image of the distribution",cmd);
	TCLAP::ValueArg<std::string> covArg("c","cov","Covariance of the distribution",true,"","covariance image of the distribution",cmd);
	TCLAP::ValueArg<std::string> resArg("o","output","Sample generated from distribution",true,"","sample output image",cmd);
    
    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default : all cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
    
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }
    
    typedef anima::SampleImageFromDistributionImageFilter <double> MainFilterType;
    
    typedef itk::ImageFileReader < MainFilterType::TInputImage > ImageReaderType;
	typedef itk::ImageFileWriter < MainFilterType::TOutputImage > ImageWriterType;

    MainFilterType::Pointer mainFilter = MainFilterType::New();
    
    ImageReaderType::Pointer meanReader = ImageReaderType::New();
    meanReader->SetFileName(meanArg.getValue());
    meanReader->Update();
    mainFilter->SetInput(0,meanReader->GetOutput());

    ImageReaderType::Pointer covReader = ImageReaderType::New();
    covReader->SetFileName(covArg.getValue());
    covReader->Update();
    mainFilter->SetInput(1,covReader->GetOutput());
    mainFilter->SetNumberOfWorkUnits(nbpArg.getValue());
    
    mainFilter->Update();
    
    std::cout << "Writing result to : " << resArg.getValue() << std::endl;
    
    ImageWriterType::Pointer outWriter = ImageWriterType::New();
    outWriter->SetInput(mainFilter->GetOutput());
    outWriter->SetFileName(resArg.getValue());
    outWriter->SetUseCompression(true);
    
    outWriter->Update();
    
    return 0;
}
