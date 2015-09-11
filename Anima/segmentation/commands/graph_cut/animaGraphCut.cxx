#include <tclap/CmdLine.h>
#include <itkNiftiImageIO.h>
#include <itkTimeProbe.h>
#include <itkImageFileReader.h>
#include "animaGraphCutFilter.h"

int main(int argc, const char** argv)
{

    // Parsing arguments
    TCLAP::CmdLine  cmd("INRIA / IRISA - VisAGeS Team", ' ',"1.0");

    TCLAP::ValueArg<std::string> inputFileT1Arg("i","input-t1","T1 input image",false,"","T1 input image",cmd);
    TCLAP::ValueArg<std::string> inputFileT2Arg("j","input-t2","T2 input image",false,"","T2 input image",cmd);
    TCLAP::ValueArg<std::string> inputFileDPArg("k","input-dp","DP input image",false,"","DP input image",cmd);
    TCLAP::ValueArg<std::string> inputFileFLAIRArg("l","input-flair","FLAIR input image",false,"","FLAIR input image",cmd);
    TCLAP::ValueArg<std::string> inputFileT1GdArg("","input-t1Gd","T1Gd input image",false,"","T1Gd input image",cmd);
    TCLAP::ValueArg<std::string> maskFileArg("m","mask","Brain mask",true,"","mask",cmd);
    TCLAP::ValueArg<std::string> sourcesMaskFileArg("","sof","Sources input image",false,"","Sources input image",cmd);
    TCLAP::ValueArg<std::string> sinksMaskFileArg("","sif","Sinks input image",false,"","Sinks input image",cmd);
    TCLAP::ValueArg<std::string> sourcesProbaFileArg("","sopf","Probabilities sources input image",false,"","Probabilities sources input image",cmd);
    TCLAP::ValueArg<std::string> sinksProbaFileArg("","sipf","Probabilities sinks input image",false,"","Probabilities sinks input image",cmd);
    TCLAP::SwitchArg verboseArg("v","verbose","verbose mode (default: false)",cmd,false);

    TCLAP::SwitchArg notUseSpecGradArg("u","no-usg","Do not use spectral gradient (default: false)",cmd,false);
    TCLAP::ValueArg<int> TLinkModeArg("","mode","Graph cut computation mode (default: 1 --> strem mode)",false,1,"graph cut computation mode",cmd);
    TCLAP::ValueArg<float> multiVarSourcesArg("","mv","Coefficient to multiply the variance value of the sources seed (default: 1)",false,1,"sources multiply variance",cmd);
    TCLAP::ValueArg<float> multiVarSinksArg("","ms","Coefficient to multiply the variance value of the seed (default: 1)",false,1,"sinks multiply variance",cmd);
    TCLAP::ValueArg<float> sigmaArg("s","sigma","sigma value (default: 0.6)",false,0.6,"sigma",cmd);
    TCLAP::ValueArg<float> alphaArg("a","alpha","Alpha value (default: 10)",false,10,"alpha",cmd);
    TCLAP::ValueArg<std::string> matrixGradArg("","mat","Spectral gradient matrix file",false,"","spectral gradient matrix file",cmd);

    TCLAP::ValueArg<std::string> outputGCFileArg("o","out","Spectral gradient matrix file",false,"","spectral gradient matrix file",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    const unsigned int Dimension = 3;

    typedef itk::Image <float,Dimension> InputImageTypeF;
    typedef InputImageTypeF::Pointer InputImagePointerF;
    typedef itk::ImageFileReader<InputImageTypeF> ReaderTypeF;

    typedef itk::Image <unsigned char,Dimension> InputImageTypeUC;
    typedef InputImageTypeUC::Pointer InputImagePointerUC;
    typedef itk::ImageFileReader<InputImageTypeUC> ReaderTypeUC;

    // Create instance of graph cut filter
    typedef anima::GraphCutFilter<InputImageTypeF,InputImageTypeUC>  FilterTypeGraphCut;
    FilterTypeGraphCut::Pointer GraphCutFilter = FilterTypeGraphCut::New();

    // Open Images
    ReaderTypeF::Pointer tmpReadF;
    ReaderTypeUC::Pointer tmpReadUC;

    if( inputFileT1Arg.getValue()!="" )
    {
        tmpReadF = ReaderTypeF::New();
        tmpReadF->SetFileName( inputFileT1Arg.getValue() );
        tmpReadF->Update();
        GraphCutFilter->SetInputImage1( tmpReadF->GetOutput() );
    }

    if( inputFileT2Arg.getValue()!="" )
    {
        tmpReadF = ReaderTypeF::New();
        tmpReadF->SetFileName( inputFileT2Arg.getValue() );
        tmpReadF->Update();
        GraphCutFilter->SetInputImage2( tmpReadF->GetOutput() );
    }

    if( inputFileDPArg.getValue()!="" )
    {
        tmpReadF = ReaderTypeF::New();
        tmpReadF->SetFileName( inputFileDPArg.getValue() );
        tmpReadF->Update();
        GraphCutFilter->SetInputImage3( tmpReadF->GetOutput() );
    }

    if( inputFileFLAIRArg.getValue()!="" )
    {
        tmpReadF = ReaderTypeF::New();
        tmpReadF->SetFileName( inputFileFLAIRArg.getValue() );
        tmpReadF->Update();
        GraphCutFilter->SetInputImage4( tmpReadF->GetOutput() );
    }

    if( inputFileT1GdArg.getValue()!="" )
    {
        tmpReadF = ReaderTypeF::New();
        tmpReadF->SetFileName( inputFileT1GdArg.getValue() );
        tmpReadF->Update();
        GraphCutFilter->SetInputImage5( tmpReadF->GetOutput() );
    }

    if( sourcesProbaFileArg.getValue()!="" )
    {
        tmpReadF = ReaderTypeF::New();
        tmpReadF->SetFileName( sourcesProbaFileArg.getValue() );
        tmpReadF->Update();
        GraphCutFilter->SetInputSeedSourcesProba( tmpReadF->GetOutput() );
    }

    if( sinksProbaFileArg.getValue()!="" )
    {
        tmpReadF = ReaderTypeF::New();
        tmpReadF->SetFileName( sinksProbaFileArg.getValue() );
        tmpReadF->Update();
        GraphCutFilter->SetInputSeedSinksProba( tmpReadF->GetOutput() );
    }


    if( maskFileArg.getValue()!="" )
    {
        tmpReadUC = ReaderTypeUC::New();
        tmpReadUC->SetFileName( maskFileArg.getValue() );
        tmpReadUC->Update();
        GraphCutFilter->SetMask( tmpReadUC->GetOutput() );
    }

    if( sourcesMaskFileArg.getValue()!="" )
    {
        tmpReadUC = ReaderTypeUC::New();
        tmpReadUC->SetFileName( sourcesMaskFileArg.getValue() );
        tmpReadUC->Update();
        GraphCutFilter->SetInputSeedSourcesMask( tmpReadUC->GetOutput() );
    }

    if( sinksMaskFileArg.getValue()!="" )
    {
        tmpReadUC = ReaderTypeUC::New();
        tmpReadUC->SetFileName( sinksMaskFileArg.getValue() );
        tmpReadUC->Update();
        GraphCutFilter->SetInputSeedSinksMask( tmpReadUC->GetOutput() );
    }

    // Set parameters
    GraphCutFilter->SetVerbose( verboseArg.getValue() );
    GraphCutFilter->SetUseSpectralGradient( !(notUseSpecGradArg.getValue()) );
    GraphCutFilter->SetTLinkMode( (TLinkMode) TLinkModeArg.getValue() );
    GraphCutFilter->SetMultiVarSources( multiVarSourcesArg.getValue() );
    GraphCutFilter->SetMultiVarSinks( multiVarSinksArg.getValue() );
    GraphCutFilter->SetAlpha( alphaArg.getValue() );
    GraphCutFilter->SetSigma( sigmaArg.getValue() );
    GraphCutFilter->SetMatrixGradFilename( matrixGradArg.getValue() );
    GraphCutFilter->SetOutputFilename( outputGCFileArg.getValue() );


    // Process
    itk::TimeProbe timer;

    timer.Start();

    try
    {
        GraphCutFilter->Update();
        GraphCutFilter->WriteOutputs();
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














