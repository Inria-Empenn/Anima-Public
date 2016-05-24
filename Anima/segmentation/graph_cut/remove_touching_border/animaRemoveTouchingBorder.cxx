#include <tclap/CmdLine.h>

#include <itkTimeProbe.h>
#include <animaReadWriteFunctions.h>
#include "animaRemoveTouchingBorderFilter.h"

int main(int argc, const char** argv)
{
    const unsigned int Dimension = 3;
    typedef itk::Image <double,Dimension> InputImageTypeD;
    typedef itk::Image <unsigned char,Dimension> InputImageTypeUC;
    typedef anima::RemoveTouchingBorderFilter<InputImageTypeUC,InputImageTypeUC,InputImageTypeUC>  FilterTypeSeg;

    // Parsing arguments
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',ANIMA_VERSION);

    // Setting up parameters

    // Input filenames
    TCLAP::ValueArg<std::string> inArg("i","inputfile","Input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","Output image",true,"","output image",cmd);
    
    TCLAP::ValueArg<std::string> maskFileArg("m","mask","Brain mask",true,"","brain mask",cmd);

    // global
    TCLAP::ValueArg<unsigned int> numThreadsArg("T","threads","Number of execution threads (default: 0 = all cores)",false,0,"number of threads",cmd);
    TCLAP::SwitchArg verboseArg("v","verbose","verbose mode (default: false)",cmd,false);
    TCLAP::ValueArg<double> tolArg("","tol","Filter tolerance (default: 0.0001)",false,0.0001,"filter tolerance",cmd);

    try
	{
	    cmd.parse(argc,argv);
	}
    catch (TCLAP::ArgException& e)
	{
	    std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
	    return(1);
	}

    FilterTypeSeg::Pointer segFilter = FilterTypeSeg::New();

    if( inArg.getValue()!="" )
	{
	    segFilter->SetInputImageSeg( anima::readImage<InputImageTypeUC>( inArg.getValue() ) );
	}

    if( outArg.getValue()!="" )
	{
	    segFilter->SetOutputNonTouchingBorderFilename(outArg.getValue());
	}
    if( maskFileArg.getValue()!="" )
	{
	    segFilter->SetMask( anima::readImage<InputImageTypeUC>( maskFileArg.getValue() ) );
	}

    // Set parameters

    segFilter->SetVerbose( verboseArg.getValue() );
    segFilter->SetNumberOfThreads( numThreadsArg.getValue() );
    segFilter->SetTol( tolArg.getValue() );

    // Process
    itk::TimeProbe timer;

    std::cout << "Remove lesions touching border..." << std::endl;
    
    timer.Start();

    try
	{
	    segFilter->Update();
	    segFilter->WriteOutputs();
	}
    catch (itk::ExceptionObject &e)
	{
	    std::cerr << e << std::endl;
	    return(1);
	}

    timer.Stop();

    std::cout << "Elapsed Time: " << timer.GetTotal()  << timer.GetUnit() << std::endl;

    return EXIT_SUCCESS;

}
