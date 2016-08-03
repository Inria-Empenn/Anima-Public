#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkRegionalMaximaImageFilter.h>

#include <tclap/CmdLine.h>

using namespace std;

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inputArg("i","inputimage","Input image",true,"","Input image",cmd);
    TCLAP::ValueArg<std::string> outputArg("o","outputimage","Output image",true,"","Output image",cmd);
    
    try
    {
	cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
	std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
	return(1);
    }

    typedef itk::Image <unsigned char,3> ImageType;
    typedef itk::Image <float,3> FloatImageType;
    typedef itk::ImageFileReader <FloatImageType> ImageReader;
    typedef itk::ImageFileWriter <ImageType> ImageWriter;

    typedef itk::RegionalMaximaImageFilter <FloatImageType,ImageType> MainFilterType;

    ImageReader::Pointer imRead = ImageReader::New();
    imRead->SetFileName(inputArg.getValue());

    try{
	imRead->Update();
    }
    catch (itk::ExceptionObject &e)
    {
	std::cerr << e << std::endl;
	return 1;
    }
    
    MainFilterType::Pointer mainFilter = MainFilterType::New();
    mainFilter->SetInput(imRead->GetOutput());
    mainFilter->SetForegroundValue(1);

    try
    {
	mainFilter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
	std::cerr << e << std::endl;
	return 1;
    }
    
    ImageWriter::Pointer tmpWriter = ImageWriter::New();
    tmpWriter->SetFileName(outputArg.getValue());
    tmpWriter->SetUseCompression(true);
    tmpWriter->SetInput(mainFilter->GetOutput());
    
    tmpWriter->Update();
	
    return 0;
}
