#include <tclap/CmdLine.h>

#include <itkTimeProbe.h>
#include <animaReadWriteFunctions.h>
#include "animaGcStremMsLesionsSegmentationFilter.h"

int main(int argc, const char** argv)
{
    const unsigned int Dimension = 3;
    typedef itk::Image <float,Dimension> InputImageTypeF;
    typedef itk::Image <unsigned char,Dimension> InputImageTypeUC;
    typedef anima::GcStremMsLesionsSegmentationFilter<InputImageTypeF>  FilterTypeSeg;

    // Parsing arguments
    TCLAP::CmdLine  cmd("INRIA / IRISA - VisAGeS Team", ' ',"1.0");

    // Setting up parameters

    // Input filenames
    TCLAP::ValueArg<std::string> inputFileT1Arg("i","input-t1","T1 input image",false,"","t1 input image",cmd);
    TCLAP::ValueArg<std::string> inputFileT2Arg("j","input-t2","T2 input image",false,"","t2 input image",cmd);
    TCLAP::ValueArg<std::string> inputFileDPArg("k","input-dp","DP input image",false,"","dp input image",cmd);
    TCLAP::ValueArg<std::string> inputFileFLAIRArg("l","input-flair","FLAIR input image",false,"","flair input image",cmd);
    TCLAP::ValueArg<std::string> inputFileT1GdArg("","input-t1Gd","T1Gd input image",false,"","t1gd input image",cmd);

    TCLAP::ValueArg<std::string> atlasFileCSFArg("x","atlas-csf","CSF atlas image",false,"","csf atlas image",cmd);
    TCLAP::ValueArg<std::string> atlasFileGMArg("y","atlas-gm","GM atlas image",false,"","gm atlas image",cmd);
    TCLAP::ValueArg<std::string> atlasFileWMArg("z","atlas-wm","WM atlas image",false,"","wm atlas image",cmd);

    TCLAP::ValueArg<std::string> maskFileArg("m","mask","Brain mask",true,"","brain mask",cmd);

    TCLAP::ValueArg<std::string> sourcesMaskFileArg("","so","Binary source input image for graph cut",false,"","sources input image",cmd);
    TCLAP::ValueArg<std::string> sinksMaskFileArg("","si","Binary sink input image for graph cut",false,"","sinks input image",cmd);

    // global
    TCLAP::ValueArg<unsigned int> numThreadsArg("T","threads","Number of execution threads (default: 0 = all cores)",false,0,"number of threads",cmd);
    TCLAP::SwitchArg verboseArg("v","verbose","verbose mode (default: false)",cmd,false);
    TCLAP::ValueArg<double> tolArg("","tol","Filter tolerance (default: 0.0001)",false,0.0001,"filter tolerance",cmd);

    TCLAP::ValueArg<unsigned int> lesionsSelectionArg("a","lesion-select-method","Lesions selection method (0: strem, 1: gcem, 2: gcem + manual graph cut, 3: manual graph cut, default: 0)",false,0,"lesions selection",cmd);

    // Automatic segmentation parameters
    TCLAP::ValueArg<unsigned int> initMethodArg("","ini","Initialisation method (0: Atlas, 1: Hierar DP, 2: Hierar FLAIR, default: 1)",false,1,"initialisation method",cmd);
    TCLAP::ValueArg<int> emIterArg("","iter","Maximum number of iterations in EM algorithm (default: 100)",false,100,"maximum number of iterations",cmd);
    TCLAP::ValueArg<double> minDistanceArg("","min-dist","Minimum distance in EM algorithm (default: 0.0001)",false,0.0001,"minimum distance",cmd);
    TCLAP::ValueArg<double> rejRatioArg("","rej","Rejection ratio for EM algorithm (default: 0.2)",false,0.2,"rejection ratio",cmd);
    TCLAP::ValueArg<int> emIter_concentrationArg("","iter-c","Maximum number of iterations between concentration step (default: 100)",false,100,"maximum number of iterations between concentration step ",cmd);
    TCLAP::SwitchArg em_before_concentrationArg("","em-c","Use EM before concentration (default: false)",cmd,false);

    TCLAP::ValueArg<float> thStremCSFArg("","th-stremCSF","Strem threshold in CSF image (default: 0.5)",false,0.5,"strem threshold CSF",cmd);
    TCLAP::ValueArg<float> thStremGMArg("","th-stremGM","Strem threshold in GM image (default: 0.5)",false,0.5,"strem threshold GM",cmd);
    TCLAP::ValueArg<float> thStremWMArg("","th-stremWM","Strem threshold in WM image (default: 0.5)",false,0.5,"strem threshold WM",cmd);
    TCLAP::ValueArg<double> fuzzyRuleMinArg("","min-fuzzy","Minimum factor in fuzzy rules (default: 1)",false,1,"Minimum in fuzzy rules ",cmd);
    TCLAP::ValueArg<double> fuzzyRuleMaxArg("","max-fuzzy","Maximum factor in fuzzy rules (default: 2)",false,2,"Maximum in fuzzy rules",cmd);

    TCLAP::ValueArg<std::string> readSolutionFileArg("","read-sol","Filename to read automatic segmentation solution",false,"","read solution filename",cmd);
    TCLAP::ValueArg<std::string> writeSolutionFileArg("","write-sol","Filename to write automatic segmentation solution",false,"","write solution filename",cmd);

    TCLAP::SwitchArg useT2Arg("","useT2","Use T2 image in automatic segmentation (default: false)",cmd,false);
    TCLAP::SwitchArg useDPArg("","useDP","Use DP image in automatic segmentation (default: false)",cmd,false);
    TCLAP::SwitchArg useFLAIRArg("","useFLAIR","Use FLAIR image in automatic segmentation(default: false)",cmd,false);

    // Graph cut parameters
    TCLAP::SwitchArg notUseSpecGradArg("","no-usg","Do not use spectral gradient (default: false)",cmd,false);
    TCLAP::ValueArg<float> multiVarSourcesArg("","mv","Coefficient to multiply the variance value of the source seeds (default: 1)",false,1,"sources multiply variance",cmd);
    TCLAP::ValueArg<float> multiVarSinksArg("","ms","Coefficient to multiply the variance value of the sink seeds (default: 1)",false,1,"sinks multiply variance",cmd);
    TCLAP::ValueArg<float> sigmaArg("","sigma","Sigma value (default: 0.6)",false,0.6,"sigma",cmd);
    TCLAP::ValueArg<float> alphaArg("","alpha","Alpha value (default: 10)",false,10,"alpha",cmd);
    TCLAP::ValueArg<std::string> matrixGradArg("","mat","Spectral gradient matrix file",false,"","spectral gradient matrix file",cmd);

    // Heuristic rules
    TCLAP::ValueArg<float> minLesionSizeArg("","ml","Minimum lesion size in mm3 (default: 0)",false,0,"minimum lesion size",cmd);
    TCLAP::SwitchArg removeBorderArg("","rb","Remove lesions that do not touch mask border (default: false)",cmd,false);
    TCLAP::ValueArg<double> intensityT2Arg("","intT2","T2 intensity (default: 0)",false,0,"T2 intensity ",cmd);
    TCLAP::ValueArg<double> intensityDPArg("","intDP","DP intensity (default: 0)",false,0,"DP intensity",cmd);
    TCLAP::ValueArg<double> intensityFLAIRArg("","intFLAIR","FLAIR intensity (default: 0)",false,0,"FLAIR intensity",cmd);
    TCLAP::ValueArg<float> thWMmapArg("","th-wm-map","Probability threshold to define white matter map (default: 0.2)",false,0.2,"threshold white matter map",cmd);
    TCLAP::ValueArg<float> contourWMRatioArg("","WMratio","Minimum percentage of WM that must surround lesions",false,0,"WM ratio contour",cmd);

    // Output filenames
    TCLAP::ValueArg<std::string> outputLesionFileArg("o","out-le","Output containing lesion segmentation",false,"","lesion segmentation output filename",cmd);
    TCLAP::ValueArg<std::string> outputCSFFileArg("","out-csf","Output containing CSF region",false,"","csf output filename",cmd);
    TCLAP::ValueArg<std::string> outputGMFileArg("","out-gm","Output containing GM region",false,"","gm output filename",cmd);
    TCLAP::ValueArg<std::string> outputWMFileArg("","out-wm","Output containing WM region",false,"","wm output filename",cmd);
    TCLAP::ValueArg<std::string> outputWholeFileArg("","out-wo","Output containing the whole brain segmentation",false,"","whole segmentation output filename",cmd);

    TCLAP::ValueArg<std::string> outputGCFileArg("","out-gc","Output containing graph cut segmentation",false,"","graph cut output filename",cmd);
    TCLAP::ValueArg<std::string> outputStremFileArg("","out-strem","Output of automatic strem computation",false,"","strem output",cmd);
    TCLAP::ValueArg<std::string> outputStremCSFFileArg("","out-strem-csf","Output of automatic csf strem computation",false,"","strem csf output",cmd);
    TCLAP::ValueArg<std::string> outputStremGMFileArg("","out-strem-gm","Output of automatic grey matter strem computation",false,"","strem gm output",cmd);
    TCLAP::ValueArg<std::string> outputStremWMFileArg("","out-strem-wm","Output of automatic white matter strem computation",false,"","strem wm output",cmd);

    TCLAP::ValueArg<std::string> outputObjectFileArg("","out-so", "Fuzzy object probability map (sources input for graph cut)",false,"","sources filename",cmd);

    TCLAP::ValueArg<std::string> outputMahaMaxiFileArg("","out-maha-maxi"," Mahalanobis maximum image output filename",false,"","mahalanobis maximum output filename",cmd);
    TCLAP::ValueArg<std::string> outputMahaMiniFileArg("","out-maha-mini","Mahalanobis minimum image output filename",false,"","mahalanobis minimum output filename",cmd);
    TCLAP::ValueArg<std::string> outputMahaCSFFileArg("","out-maha-csf","Mahalanobis CSF image output filename",false,"","mahalanobis csf output filename",cmd);
    TCLAP::ValueArg<std::string> outputMahaGMFileArg("","out-maha-gm","Mahalanobis GM image output filename",false,"","mahalanobis gm output filename",cmd);
    TCLAP::ValueArg<std::string> outputMahaWMFileArg("","out-maha-wm","Mahalanobis WM image output filename",false,"","mahalanobis wm output filename",cmd);

    TCLAP::ValueArg<std::string> outputHyperIntensity1FileArg("","out-hyper-im1","Hyper-intensity map in image 1 filename",false,"","output hyper-intensity filename 1",cmd);
    TCLAP::ValueArg<std::string> outputHyperIntensity2FileArg("","out-hyper-im2","Hyper-intensity map in image 2 filename",false,"","output hyper-intensity filename 2",cmd);

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

    if( inputFileT1Arg.getValue()!="" )
    {
        segFilter->SetInputImageT1( anima::readImage<InputImageTypeF>( inputFileT1Arg.getValue() ) );
    }

    if( inputFileT2Arg.getValue()!="" )
    {
        segFilter->SetInputImageT2( anima::readImage<InputImageTypeF>( inputFileT2Arg.getValue() ) );
    }

    if( inputFileDPArg.getValue()!="" )
    {
        segFilter->SetInputImageDP( anima::readImage<InputImageTypeF>( inputFileDPArg.getValue() ) );
    }

    if( inputFileFLAIRArg.getValue()!="" )
    {
        segFilter->SetInputImageFLAIR( anima::readImage<InputImageTypeF>( inputFileFLAIRArg.getValue() ) );
    }

    if( inputFileT1GdArg.getValue()!="" )
    {
        segFilter->SetInputImageT1Gd( anima::readImage<InputImageTypeF>( inputFileT1GdArg.getValue() ) );
    }

    if( atlasFileCSFArg.getValue()!="" )
    {
        segFilter->SetInputCSFAtlas( anima::readImage<InputImageTypeF>( atlasFileCSFArg.getValue() ) );
    }

    if( atlasFileGMArg.getValue()!="" )
    {
        segFilter->SetInputGMAtlas( anima::readImage<InputImageTypeF>( atlasFileGMArg.getValue() ) );
    }

    if( atlasFileWMArg.getValue()!="" )
    {
        segFilter->SetInputWMAtlas( anima::readImage<InputImageTypeF>( atlasFileWMArg.getValue() ) );
    }

    if( maskFileArg.getValue()!="" )
    {
        segFilter->SetMask( anima::readImage<InputImageTypeF>( maskFileArg.getValue() ) );
    }

    if( sourcesMaskFileArg.getValue()!="" )
    {
        segFilter->SetSourcesMask( anima::readImage<InputImageTypeUC>( sourcesMaskFileArg.getValue()) );
    }

    if( sinksMaskFileArg.getValue()!="" )
    {
        segFilter->SetSinksMask( anima::readImage<InputImageTypeUC>( sinksMaskFileArg.getValue()) );
    }

    // Set parameters

    segFilter->SetLesionSegmentationType( (LesionSegmentationType) lesionsSelectionArg.getValue() );
    segFilter->SetVerbose( verboseArg.getValue() );
    segFilter->SetNumberOfThreads( numThreadsArg.getValue() );
    segFilter->SetTol( tolArg.getValue() );

    segFilter->SetUseT2( useT2Arg.getValue() );
    segFilter->SetUseDP( useDPArg.getValue() );
    segFilter->SetUseFLAIR( useFLAIRArg.getValue() );
    segFilter->SetInitMethodType( (InitializationType) initMethodArg.getValue() );
    segFilter->SetEmIter( emIterArg.getValue() );
    segFilter->SetMinDistance( minDistanceArg.getValue() );
    segFilter->SetRejRatio( rejRatioArg.getValue());
    segFilter->SetEmIter_concentration( emIter_concentrationArg.getValue() );
    segFilter->SetEM_before_concentration( em_before_concentrationArg.getValue() );
    segFilter->SetMahalanobisThCSF( thStremCSFArg.getValue() );
    segFilter->SetMahalanobisThGM( thStremGMArg.getValue() );
    segFilter->SetMahalanobisThWM( thStremWMArg.getValue() );
    segFilter->SetFuzzyRuleMin( fuzzyRuleMinArg.getValue() );
    segFilter->SetFuzzyRuleMax( fuzzyRuleMaxArg.getValue() );
    segFilter->SetSolutionReadFilename( readSolutionFileArg.getValue() );
    segFilter->SetSolutionWriteFilename( writeSolutionFileArg.getValue() );

    segFilter->SetUseSpecGrad( !(notUseSpecGradArg.getValue()) );
    segFilter->SetMultiVarSources( multiVarSourcesArg.getValue() );
    segFilter->SetMultiVarSinks( multiVarSinksArg.getValue() );
    segFilter->SetAlpha( alphaArg.getValue() );
    segFilter->SetSigma( sigmaArg.getValue() );
    segFilter->SetMatrixGradFilename( matrixGradArg.getValue() );

    segFilter->SetIntensityT2Factor( intensityT2Arg.getValue() );
    segFilter->SetIntensityDPFactor( intensityDPArg.getValue() );
    segFilter->SetIntensityFLAIRFactor( intensityFLAIRArg.getValue() );
    segFilter->SetRemoveBorder( removeBorderArg.getValue() );
    segFilter->SetMinLesionsSize( minLesionSizeArg.getValue() );
    segFilter->SetThresoldWMmap( thWMmapArg.getValue());
    segFilter->SetRatioContourWM( contourWMRatioArg.getValue());

    segFilter->SetOutputLesionFilename( outputLesionFileArg.getValue() );
    segFilter->SetOutputCSFFilename( outputCSFFileArg.getValue() );
    segFilter->SetOutputGMFilename( outputGMFileArg.getValue() );
    segFilter->SetOutputWMFilename( outputWMFileArg.getValue() );
    segFilter->SetOutputWholeFilename( outputWholeFileArg.getValue() );

    segFilter->SetOutputGCFilename( outputGCFileArg.getValue() );
    segFilter->SetOutputStremFilename( outputStremFileArg.getValue() );
    segFilter->SetOutputStremCSFFilename( outputStremCSFFileArg.getValue() );
    segFilter->SetOutputStremGMFilename( outputStremGMFileArg.getValue() );
    segFilter->SetOutputStremWMFilename( outputStremWMFileArg.getValue() );

    segFilter->SetOutputFuzzyObjectFilename( outputObjectFileArg.getValue() );
    segFilter->SetOutputMahaMaximumFilename( outputMahaMaxiFileArg.getValue() );
    segFilter->SetOutputMahaMinimumFilename( outputMahaMiniFileArg.getValue() );
    segFilter->SetOutputMahaCSFFilename( outputMahaCSFFileArg.getValue() );
    segFilter->SetOutputMahaGMFilename( outputMahaGMFileArg.getValue() );
    segFilter->SetOutputMahaWMFilename( outputMahaWMFileArg.getValue() );

    segFilter->SetOutputIntensityImage1Filename( outputHyperIntensity1FileArg.getValue() );
    segFilter->SetOutputIntensityImage2Filename( outputHyperIntensity2FileArg.getValue() );

    // Process
    itk::TimeProbe timer;

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
