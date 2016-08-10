#include <tclap/CmdLine.h>

#include <itkTimeProbe.h>
#include <animaReadWriteFunctions.h>
#include "animaRemoveTouchingBorderFilter.h"

int main(int argc, const char** argv)
{
    const unsigned int Dimension = 3;
    typedef itk::Image <unsigned char,Dimension> MaskImageType;
    typedef itk::Image <unsigned int,Dimension> InputImageType;
    typedef anima::RemoveTouchingBorderFilter<InputImageType,MaskImageType,MaskImageType>  FilterTypeSeg;
    
    // Parsing arguments
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',ANIMA_VERSION);

    // Setting up parameters

    // Input filenames
    TCLAP::ValueArg<std::string> inArg("i","inputfile","Input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","Output image",true,"","output image",cmd);
    
    TCLAP::ValueArg<std::string> maskFileArg("m","mask","Brain mask",true,"","brain mask",cmd);

    TCLAP::SwitchArg noContourDetectionArg("c","noContour","Avoid contour detection (default: false)",cmd,false);
    TCLAP::SwitchArg rmNonTouchingArg("","rmNonTouching","Remove non touching components instead of touching components (default: false)",cmd,false);
    TCLAP::SwitchArg labeledImageArg("l","labeledImage","Input image is a labeled image (default: false)",cmd,false);

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
    
    if(labeledImageArg.getValue())
    {
	typedef anima::RemoveTouchingBorderFilter<MaskImageType,MaskImageType,MaskImageType>  FilterTypeSeg;

	FilterTypeSeg::Pointer segFilter = FilterTypeSeg::New();
    
	if( inArg.getValue()!="" )
	{
	    segFilter->SetInputImageSeg( anima::readImage<MaskImageType>( inArg.getValue() ) );
	}
    }
    else
    {
	if( inArg.getValue()!="" )
	{
	    segFilter->SetInputImageSeg( anima::readImage<InputImageType>( inArg.getValue() ) );
	}
    }
	
    if( outArg.getValue()!="" )
    {
        if( rmNonTouchingArg.getValue() )
        {
            segFilter->SetOutputTouchingBorderFilename(outArg.getValue());
        }
        else
        {
            segFilter->SetOutputNonTouchingBorderFilename(outArg.getValue());
        }
    }
    if( maskFileArg.getValue()!="" )
    {
        segFilter->SetMask( anima::readImage<MaskImageType>( maskFileArg.getValue() ) );
    }

    // Set parameters

    segFilter->SetLabeledImage( labeledImageArg.getValue() );
    segFilter->SetNoContour( noContourDetectionArg.getValue() );
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
