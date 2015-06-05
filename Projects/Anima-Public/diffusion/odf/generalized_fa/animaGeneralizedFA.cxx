#include <animaGeneralizedFAImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
  TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',"1.0");
	
  TCLAP::ValueArg<std::string> inArg("i","inputodf","ODF volume",true,"","ODF volume",cmd);
	TCLAP::ValueArg<std::string> resArg("o","output","Result image",true,"","result GFA image",cmd);
	
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
	
    typedef anima::GeneralizedFAImageFilter <float> MainFilterType;
	typedef itk::ImageFileReader < MainFilterType::TInputImage > ImageReaderType;
	typedef itk::ImageFileWriter < MainFilterType::TOutputImage > ImageWriterType;
	
	MainFilterType::Pointer mainFilter = MainFilterType::New();
	
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
	imageReader->SetFileName(inArg.getValue());
	imageReader->Update();
	
	mainFilter->SetInput(imageReader->GetOutput());
	mainFilter->SetNumberOfThreads(nbpArg.getValue());
	
	mainFilter->Update();
	
	ImageWriterType::Pointer imageWriter = ImageWriterType::New();
	imageWriter->SetInput(mainFilter->GetOutput());
	imageWriter->SetFileName(resArg.getValue());
	imageWriter->SetUseCompression(true);
	
	imageWriter->Update();
	
	return 0;
}
