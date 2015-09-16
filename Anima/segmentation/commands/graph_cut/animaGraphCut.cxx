#include <tclap/CmdLine.h>

#include <itkTimeProbe.h>
#include <animaReadWriteFunctions.h>
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
    TCLAP::ValueArg<double> multiVarSourcesArg("","mv","Coefficient to multiply the variance value of the sources seed (default: 1)",false,1,"sources multiply variance",cmd);
    TCLAP::ValueArg<double> multiVarSinksArg("","ms","Coefficient to multiply the variance value of the seed (default: 1)",false,1,"sinks multiply variance",cmd);
    TCLAP::ValueArg<double> sigmaArg("s","sigma","sigma value (default: 0.6)",false,0.6,"sigma",cmd);
    TCLAP::ValueArg<double> alphaArg("a","alpha","Alpha value (default: 10)",false,10,"alpha",cmd);
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
    typedef itk::Image <double,Dimension> InputImageTypeD;
    typedef itk::Image <unsigned char,Dimension> InputImageTypeUC;

    // Create instance of graph cut filter
    typedef anima::GraphCutFilter<InputImageTypeD,InputImageTypeUC>  FilterTypeGraphCut;
    FilterTypeGraphCut::Pointer GraphCutFilter = FilterTypeGraphCut::New();

    if( inputFileT1Arg.getValue()!="" )
    {
        GraphCutFilter->SetInputImage1( anima::readImage<InputImageTypeD>(inputFileT1Arg.getValue()) );
    }

    if( inputFileT2Arg.getValue()!="" )
    {
        GraphCutFilter->SetInputImage2( anima::readImage<InputImageTypeD>(inputFileT2Arg.getValue()) );
    }

    if( inputFileDPArg.getValue()!="" )
    {
        GraphCutFilter->SetInputImage3( anima::readImage<InputImageTypeD>(inputFileDPArg.getValue()) );
    }

    if( inputFileFLAIRArg.getValue()!="" )
    {
        GraphCutFilter->SetInputImage4( anima::readImage<InputImageTypeD>(inputFileFLAIRArg.getValue()) );
    }

    if( inputFileT1GdArg.getValue()!="" )
    {
        GraphCutFilter->SetInputImage5( anima::readImage<InputImageTypeD>(inputFileT1GdArg.getValue()) );
    }

    if( sourcesProbaFileArg.getValue()!="" )
    {
        GraphCutFilter->SetInputSeedSourcesProba( anima::readImage<InputImageTypeD>(sourcesProbaFileArg.getValue()) );
    }

    if( sinksProbaFileArg.getValue()!="" )
    {
        GraphCutFilter->SetInputSeedSinksProba( anima::readImage<InputImageTypeD>(sinksProbaFileArg.getValue()) );
    }

    if( maskFileArg.getValue()!="" )
    {
        GraphCutFilter->SetMask( anima::readImage<InputImageTypeUC>(maskFileArg.getValue()) );
    }

    if( sourcesMaskFileArg.getValue()!="" )
    {
        GraphCutFilter->SetInputSeedSourcesMask( anima::readImage<InputImageTypeUC>(sourcesMaskFileArg.getValue()) );
    }

    if( sinksMaskFileArg.getValue()!="" )
    {
        GraphCutFilter->SetInputSeedSinksMask( anima::readImage<InputImageTypeUC>(sinksMaskFileArg.getValue()) );
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














