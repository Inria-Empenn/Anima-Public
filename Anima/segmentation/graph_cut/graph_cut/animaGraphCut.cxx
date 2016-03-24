#include <tclap/CmdLine.h>

#include <itkTimeProbe.h>
#include <animaReadWriteFunctions.h>
#include "animaGraphCutFilter.h"

int main(int argc, const char** argv)
{

    // Parsing arguments
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',ANIMA_VERSION);

    TCLAP::MultiArg<std::string> inputArg("i","input","Input images to segment (up to 5, e.g. T1, T2, DP, FLAIR, T1-Gd)",
                                          true,"input images",cmd);

    TCLAP::ValueArg<std::string> maskFileArg("m","mask","Brain mask",false,"","brain mask",cmd);
    TCLAP::ValueArg<std::string> sourcesMaskFileArg("","sof","Sources input image",true,"","Sources input image",cmd);
    TCLAP::ValueArg<std::string> sinksMaskFileArg("","sif","Sinks input image",true,"","Sinks input image",cmd);
    TCLAP::ValueArg<std::string> sourcesProbaFileArg("","sopf","Probabilities sources input image",false,"","Probabilities sources input image",cmd);
    TCLAP::ValueArg<std::string> sinksProbaFileArg("","sipf","Probabilities sinks input image",false,"","Probabilities sinks input image",cmd);
    TCLAP::SwitchArg verboseArg("V","verbose","verbose mode (default: false)",cmd,false);

    TCLAP::SwitchArg notUseSpecGradArg("U","no-usg","Do not use spectral gradient (default: false)",cmd,false);
    TCLAP::ValueArg<int> TLinkModeArg("","mode","Graph cut computation mode (0: single Gaussian, 1: STREM mode, default: 0)",false,0,"graph cut computation mode",cmd);
    TCLAP::ValueArg<double> multiVarSourcesArg("","mv","Coefficient to multiply the variance value of the sources seed (default: 1)",false,1,"sources multiply variance",cmd);
    TCLAP::ValueArg<double> multiVarSinksArg("","ms","Coefficient to multiply the variance value of the seed (default: 1)",false,1,"sinks multiply variance",cmd);
    TCLAP::ValueArg<double> sigmaArg("s","sigma","sigma value (default: 0.6)",false,0.6,"sigma",cmd);
    TCLAP::ValueArg<double> alphaArg("a","alpha","Alpha value (default: 10)",false,10,"alpha",cmd);
    TCLAP::ValueArg<std::string> matrixGradArg("","mat","Spectral gradient matrix file",false,"","spectral gradient matrix file",cmd);

    TCLAP::ValueArg<std::string> outputGCFileArg("o","out","Output segmentation",true,"","output segmentation",cmd);

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

    for (unsigned int i = 0;i < inputArg.getValue().size();++i)
        GraphCutFilter->SetInputImage(i,anima::readImage<InputImageTypeD>(inputArg.getValue()[i]));

    if( sourcesProbaFileArg.getValue()!="" )
    {
        GraphCutFilter->SetInputSeedSourcesProba( anima::readImage<InputImageTypeD>(sourcesProbaFileArg.getValue()) );
    }

    if( sinksProbaFileArg.getValue()!="" )
    {
        GraphCutFilter->SetInputSeedSinksProba( anima::readImage<InputImageTypeD>(sinksProbaFileArg.getValue()) );
    }

    if( maskFileArg.getValue()!="" )
        GraphCutFilter->SetMask( anima::readImage<InputImageTypeUC>(maskFileArg.getValue()) );
    else
    {
        InputImageTypeUC::Pointer maskImage = InputImageTypeUC::New();
        maskImage->SetRegions(GraphCutFilter->GetInput(0)->GetLargestPossibleRegion());
        maskImage->SetSpacing(GraphCutFilter->GetInput(0)->GetSpacing());
        maskImage->SetOrigin(GraphCutFilter->GetInput(0)->GetOrigin());
        maskImage->SetDirection(GraphCutFilter->GetInput(0)->GetDirection());

        maskImage->Allocate();
        maskImage->FillBuffer(1);

        GraphCutFilter->SetMask(maskImage);
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
