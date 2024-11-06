#include <tclap/CmdLine.h>

#include <animaReadWriteFunctions.h>
#include <itkTimeProbe.h>
#include "animaConnectedComponentsMetricsFilter.h"


struct arguments
{
    std::string ref, test, outputCCref, outputCCtest, outputDiffVoxel, outputDiffRef, outputDiffTest,
            outputEvolutionRef, outputEvolutionTest, outputTextVolRef, outputTextVolTest, outputTextDiff;
    unsigned int pthread;
    double min, alpha, beta;
    bool full, verbose, detection;
};

template <unsigned int Dimension>
void
connectedComponent(const arguments &args)
{
    typedef itk::Image <unsigned short,Dimension> ImageType;

    // Create instance of metrics 3D filter
    typedef anima::ConnectedComponentsMetricsFilter<ImageType,ImageType>  ConnectedComponentsMetricsFilterType;
    typename ConnectedComponentsMetricsFilterType::Pointer ConnectedComponentsMetricsFilter = ConnectedComponentsMetricsFilterType::New();

    // Open Images
    ConnectedComponentsMetricsFilter->SetInputReference( anima::readImage<ImageType>( args.ref ) );
    ConnectedComponentsMetricsFilter->SetInputTest( anima::readImage<ImageType>( args.test ) );

    ConnectedComponentsMetricsFilter->SetOverlapDetectionType( args.detection );
    ConnectedComponentsMetricsFilter->SetAlpha( args.alpha );
    ConnectedComponentsMetricsFilter->SetBeta( args.beta );
    ConnectedComponentsMetricsFilter->SetMinSizeMM3( args.min );
    ConnectedComponentsMetricsFilter->SetFullyConnected( args.full );
    ConnectedComponentsMetricsFilter->SetVerbose( args.verbose );
    ConnectedComponentsMetricsFilter->SetDimension( Dimension );
    if(args.pthread > 0)
        ConnectedComponentsMetricsFilter->SetNumberOfWorkUnits(args.pthread);

    ConnectedComponentsMetricsFilter->SetOutputFilenameCCRef( args.outputCCref );
    ConnectedComponentsMetricsFilter->SetOutputFilenameCCTest( args.outputCCtest );

    ConnectedComponentsMetricsFilter->SetOutputFilenameDiffVoxelWise( args.outputDiffVoxel );
    ConnectedComponentsMetricsFilter->SetOutputFilenameDiffTest( args.outputDiffTest );
    ConnectedComponentsMetricsFilter->SetOutputFilenameDiffRef( args.outputDiffRef );

    ConnectedComponentsMetricsFilter->SetOutputFilenameTextVolTest( args.outputTextVolTest );
    ConnectedComponentsMetricsFilter->SetOutputFilenameTextVolRef( args.outputTextVolRef );
    ConnectedComponentsMetricsFilter->SetOutputFilenameTextMetrics( args.outputTextDiff );

    ConnectedComponentsMetricsFilter->SetOutputFilenameEvolutionTest( args.outputEvolutionTest );
    ConnectedComponentsMetricsFilter->SetOutputFilenameEvolutionRef( args.outputEvolutionRef );

    // Process
    itk::TimeProbe timer;

    timer.Start();

    ConnectedComponentsMetricsFilter->Update();
    ConnectedComponentsMetricsFilter->WriteOutputs();
    timer.Stop();

    std::cout << "Elapsed Time: " << timer.GetTotal()  << timer.GetUnit() << std::endl;
}

void
retrieveDimension(const arguments &args, itk::ImageIOBase::Pointer imageIO)
{
    unsigned int nbDim = imageIO->GetNumberOfDimensions();

    switch(nbDim)
    {
    case 2:
        connectedComponent<2>(args);
        break;
    case 3:
        connectedComponent<3>(args);
        break;
    case 4:
        connectedComponent<4>(args);
        break;
    default:
        itk::ExceptionObject excp(__FILE__, __LINE__, "Number of dimension not supported.", ITK_LOCATION);
        throw excp;
    }
}

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',"1.0");

    TCLAP::ValueArg<std::string> refArg("r","ref","Reference input image",true,"","input reference image",cmd);
    TCLAP::ValueArg<std::string> testArg("t","test","Test input image",true,"","input test image",cmd);

    TCLAP::ValueArg<std::string> outCCRefArg("","output-cc-ref","Output connected components for reference image",false,"","reference connected components",cmd);
    TCLAP::ValueArg<std::string> outCCTestArg("","output-cc-test","Output connected components for test image",false,"","reference connected components",cmd);    
    TCLAP::ValueArg<std::string> outDiffVoxelArg("","output-diff-voxel","Output difference image voxel wise",false,"","voxel wise diff",cmd);
    TCLAP::ValueArg<std::string> outDiffRefArg("","output-diff-ref","Output difference object wise for reference image",false,"","object wise diff",cmd);
    TCLAP::ValueArg<std::string> outDiffTestArg("","output-diff-test","Output difference object wise for test image",false,"","object wise diff",cmd);
    TCLAP::ValueArg<std::string> outEvolutionRefArg("","output-evol-ref","Output evolution for reference image",false,"","evolution reference",cmd);
    TCLAP::ValueArg<std::string> outEvolutionTestArg("","output-evol-test","Output evolution for test image",false,"","evolution test",cmd);

    TCLAP::ValueArg<std::string> outTextVolRefArg("","output-vol-ref","Output text file containing reference image volume information",false,"","volume reference information",cmd);
    TCLAP::ValueArg<std::string> outTextVolTestArg("","output-vol-test","Output text file containing test image volume information",false,"","volume test information",cmd);
    TCLAP::ValueArg<std::string> outTextDiffArg("","output-metrics","Output text file containing metrics results",false,"","metrics information",cmd);

    TCLAP::ValueArg<double> alphaArg("a","alpha","overlap factor",false,0,"overlap factor",cmd);
    TCLAP::ValueArg<double> betaArg("b","beta","size evolution factor",false,0,"size evolution factor",cmd);
    TCLAP::ValueArg<double> minSizeArg("m","minsize","minimal component size in mm3",false,0,"minimal component size",cmd);
    TCLAP::SwitchArg fullConnectArg("F","full-connect","Use 26-connectivity instead of 6-connectivity",cmd,false);
    TCLAP::SwitchArg overlapDetectionArg("d","overlap-detect","use overlap detection type",cmd,false);
    TCLAP::SwitchArg verboseArg("v","verbose","print information",cmd,false);
    TCLAP::ValueArg<unsigned int> numThreadsArg("T","threads","Number of execution threads (default: 0 = all cores)",false,0,"number of threads",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    // Find out the type of the image in reference file
    itk::ImageIOBase::Pointer imageIOref = itk::ImageIOFactory::CreateImageIO(refArg.getValue().c_str(),itk::IOFileModeEnum::ReadMode);
    if( !imageIOref )
    {
        std::cerr << "Itk could not find suitable IO factory for the input" << std::endl;
        return EXIT_FAILURE;
    }

    // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
    imageIOref->SetFileName(refArg.getValue());
    imageIOref->ReadImageInformation();


    // Find out the type of the image in test file
    itk::ImageIOBase::Pointer imageIOtest = itk::ImageIOFactory::CreateImageIO(testArg.getValue().c_str(),itk::IOFileModeEnum::ReadMode);
    if( !imageIOtest )
    {
        std::cerr << "Itk could not find suitable IO factory for the input" << std::endl;
        return EXIT_FAILURE;
    }

    // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
    imageIOtest->SetFileName(testArg.getValue());
    imageIOtest->ReadImageInformation();

    if(imageIOref->GetNumberOfDimensions() != imageIOtest->GetNumberOfDimensions())
    {
        std::cerr << "Input images do not have the same sizes" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout<<"\npreparing filter...\n";

    arguments args;

    args.ref = refArg.getValue();
    args.test = testArg.getValue();

    args.outputCCref = outCCRefArg.getValue();
    args.outputCCtest = outCCTestArg.getValue();

    args.outputDiffVoxel = outDiffVoxelArg.getValue();
    args.outputDiffRef = outDiffRefArg.getValue();
    args.outputDiffTest = outDiffTestArg.getValue();

    args.outputEvolutionRef = outEvolutionRefArg.getValue();
    args.outputEvolutionTest = outEvolutionTestArg.getValue();

    args.outputTextVolRef = outTextVolRefArg.getValue();
    args.outputTextVolTest = outTextVolTestArg.getValue();
    args.outputTextDiff = outTextDiffArg.getValue();

    args.detection = overlapDetectionArg.getValue();
    args.alpha = alphaArg.getValue();
    args.beta = betaArg.getValue();
    args.verbose = verboseArg.getValue();

    args.min = minSizeArg.getValue();
    args.full= fullConnectArg.getValue();
    args.pthread = numThreadsArg.getValue();

    try
    {
        retrieveDimension(args, imageIOref);
    }
    catch ( itk::ExceptionObject & err )
    {
        std::cerr << "Itk cannot concatenate, be sure to use valid arguments..." << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }
}

