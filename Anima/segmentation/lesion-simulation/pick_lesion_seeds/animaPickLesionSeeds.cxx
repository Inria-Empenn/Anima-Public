#include <animaPickLesionSeedImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inArg("i","inputfile","Input probability image",true,"","input probability image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","Output seeds image",true,"","output seeds image",cmd);
	
    TCLAP::ValueArg<double> minDistArg("d","dist","Minimal distance between seeds (default: 10 mm)",false,10,"minimal distance between seeds",cmd);
    TCLAP::ValueArg<unsigned short> numSeedsArg("n","numseeds","Number of seeds (default: 1)",false,1,"number of seeds",cmd);
	
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }
	
    typedef anima::PickLesionSeedImageFilter <itk::Image <double, 3>, itk::Image <unsigned short, 3> > MainFilterType;
    typedef MainFilterType::InputImageType InputImageType;
    typedef MainFilterType::OutputImageType OutputImageType;
	typedef itk::ImageFileReader <InputImageType> itkImageReader;
	typedef itk::ImageFileWriter <OutputImageType> itkImageWriter;
	
	itkImageReader::Pointer imRead = itkImageReader::New();
	imRead->SetFileName(inArg.getValue());
	imRead->Update();
    
    MainFilterType::Pointer mainFilter = MainFilterType::New();
    mainFilter->SetInput(imRead->GetOutput());
    
    mainFilter->SetNumberOfSeeds(numSeedsArg.getValue());
    mainFilter->SetProximityThreshold(minDistArg.getValue());
    
    try
    {
        mainFilter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return 1;
    }
    
	itkImageWriter::Pointer tmpWriter = itkImageWriter::New();
	tmpWriter->SetFileName(outArg.getValue());
	tmpWriter->SetUseCompression(true);
	tmpWriter->SetInput(mainFilter->GetOutput());
	
	tmpWriter->Update();
	
	return 0;
}
