#include <tclap/CmdLine.h>

#include <itkTimeProbe.h>
#include <animaReadWriteFunctions.h>
#include "animaRemoveTouchingBorderFilter.h"

int main(int argc, const char** argv)
{
    const unsigned int Dimension = 3;
    typedef itk::Image <unsigned char,Dimension> InputImageType;
    typedef itk::ImageFileReader <InputImageType> ReaderType;
    typedef itk::ImageRegionConstIterator <InputImageType> InputImageIterator;
    
    // Parsing arguments
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',ANIMA_VERSION);

    // Setting up parameters

    // Input filenames
    TCLAP::ValueArg<std::string> inArg("i","inputfile","Input image",true,"","input image",cmd);
    
    try
    {
	cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
	std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
	return(1);
    }

    ReaderType::Pointer inputImageReader = ReaderType::New();
    inputImageReader->SetFileName(inputArg.getValue());

    try
    {
	inputImageReader->Update();
    }
    catch (itk::ExceptionObject &e)
    {
	std::cerr << e << std::endl;
	return 1;
    }
    
    InputImageIterator inputIt(inputImageReader->GetOutput(),inputImageReader->GetOutput()->GetLargestPossibleRegion());

    unsigned int cpt=0;
    while(!inputIt.IsAtEnd())
    {
	if( inputIt.Get()!=0 )
	{
	    ++cpt;
	}
	++inputIt;
    }

    std::cout << cpt << std::endl;
    
    return EXIT_SUCCESS;

}
