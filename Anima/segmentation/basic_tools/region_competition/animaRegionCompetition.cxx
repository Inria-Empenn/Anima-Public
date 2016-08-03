#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkRegionCompetitionImageFilter.h>

#include <tclap/CmdLine.h>

using namespace std;

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inputArg("i","inputimage","Input image",true,"","Input image",cmd);
    TCLAP::ValueArg<std::string> labeledArg("l", "labeledImage", "Labeled Image", true, "", "Labeled Image", cmd);
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

    typedef itk::Image <unsigned char,3> InputImageType;
    typedef itk::Image <unsigned short, 3> OutputImageType; 
    typedef itk::ImageFileReader <InputImageType> InputImageReader;
    typedef itk::ImageFileReader <OutputImageType> LabeledImageReader;
    typedef itk::ImageFileWriter <OutputImageType> ImageWriter;

    typedef itk::RegionCompetitionImageFilter <InputImageType,OutputImageType> MainFilterType;

    InputImageReader::Pointer inputRead = InputImageReader::New();
    inputRead->SetFileName(inputArg.getValue());

    try{
	inputRead->Update();
    }
    catch (itk::ExceptionObject &e)
    {
	std::cerr << e << std::endl;
	return 1;
    }

    LabeledImageReader::Pointer labeledRead = LabeledImageReader::New();
    labeledRead->SetFileName(labeledArg.getValue());

    try{
	labeledRead->Update();
    }
    catch (itk::ExceptionObject &e)
    {
	std::cerr << e << std::endl;
	return 1;
    }
    
    MainFilterType::Pointer mainFilter = MainFilterType::New();
    mainFilter->SetInput(inputRead->GetOutput());
    mainFilter->SetInputLabels(labeledRead->GetOutput());

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
