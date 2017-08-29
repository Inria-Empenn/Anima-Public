#include <tclap/CmdLine.h>

#include <animaPyramidalBlockMatchingBridge.h>

#include <itkTimeProbe.h>
#include <itkImageFileReader.h>
#include <itkTransformFileReader.h>

int main(int argc, const char** argv)
{
    const unsigned int Dimension = 3;

    typedef itk::Image <float,Dimension> InputImageType;
    typedef anima::PyramidalBlockMatchingBridge <Dimension> PyramidBMType;
    typedef anima::BaseTransformAgregator < Dimension > AgregatorType;
    typedef itk::AffineTransform<AgregatorType::ScalarType,Dimension> AffineTransformType;
    typedef AffineTransformType::Pointer AffineTransformPointer;

    PyramidBMType::Pointer matcher = PyramidBMType::New();

    // Parsing arguments
    TCLAP::CmdLine  cmd("INRIA / IRISA - VisAGeS Team", ' ',ANIMA_VERSION);

    // Setting up parameters
    TCLAP::ValueArg<std::string> fixedArg("r","refimage","Fixed image",true,"","fixed image",cmd);
    TCLAP::ValueArg<std::string> movingArg("m","movingimage","Moving image",true,"","moving image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputimage","Output (registered) image",true,"","output image",cmd);
    TCLAP::ValueArg<unsigned int> outTrTypeArg("","ot","Output transformation type (0: rigid, 1: translation, 2: affine, 3: anisotropic_sim, default: 0)",false,0,"output transformation type",cmd);

    TCLAP::ValueArg<std::string> initialTransformArg("i","inittransform","Initial transformation",false,"","initial transform",cmd);
    TCLAP::ValueArg<std::string> outputTransformArg("O","outtransform","Output transformation",false,"","output transform",cmd);

    TCLAP::ValueArg<std::string> blockMaskArg("M","mask-im","Mask image for block generation",false,"","block mask image",cmd);
    TCLAP::ValueArg<unsigned int> blockSizeArg("","bs","Block size (default: 5)",false,5,"block size",cmd);
    TCLAP::ValueArg<unsigned int> blockSpacingArg("","sp","Block spacing (default: 5)",false,5,"block spacing",cmd);
    TCLAP::ValueArg<float> stdevThresholdArg("s","stdev","Threshold block standard deviation (default: 5)",false,5,"block minimal standard deviation",cmd);
    TCLAP::ValueArg<double> percentageKeptArg("k","per-kept","Percentage of blocks with the highest variance kept (default: 0.8)",false,0.8,"percentage of blocks kept",cmd);

    TCLAP::ValueArg<unsigned int> blockTransfoArg("t","in-transform","Transformation computed between blocks (0: translation, 1: rigid, 2: affine, 3: directional affine, default: 0)",false,0,"transformation between blocks",cmd);
    TCLAP::ValueArg<unsigned int> directionArg("d","dir","Affine direction for directional transform output (default: 1 = Y axis)",false,1,"direction of directional affine",cmd);
    TCLAP::ValueArg<unsigned int> blockMetricArg("","metric","Similarity metric between blocks (0: squared correlation coefficient, 1: correlation coefficient, 2: mean squares, default: 0)",false,0,"similarity metric",cmd);
    TCLAP::ValueArg<unsigned int> optimizerArg("","opt","Optimizer for optimal block search (0: Exhaustive, 1: Bobyqa, default: 1)",false,1,"optimizer",cmd);

    TCLAP::ValueArg<unsigned int> maxIterationsArg("","mi","Maximum block match iterations (default: 10)",false,10,"maximum iterations",cmd);
    TCLAP::ValueArg<float> minErrorArg("","me","Minimal distance between consecutive estimated transforms (default: 0.01)",false,0.01,"minimal distance between transforms",cmd);

    TCLAP::ValueArg<unsigned int> optimizerMaxIterationsArg("","oi","Maximum iterations for local optimizer (default: 100)",false,100,"maximum local optimizer iterations",cmd);
    TCLAP::SwitchArg initIdentityArg("I","init-identity", "If no input transformation is given, initialize to identity (otherwise uses center of mass alignment)", cmd, false);

    TCLAP::ValueArg<double> searchRadiusArg("","sr","Search radius in pixels (exhaustive search window, rho start for bobyqa, default: 2)",false,2,"optimizer search radius",cmd);
    TCLAP::ValueArg<double> searchAngleRadiusArg("","sar","Search angle radius in degrees (rho start for bobyqa, default: 5)",false,5,"optimizer search angle radius",cmd);
    TCLAP::ValueArg<double> searchScaleRadiusArg("","scr","Search scale radius (rho start for bobyqa, default: 0.1)",false,0.1,"optimizer search scale radius",cmd);
    TCLAP::ValueArg<double> finalRadiusArg("","fr","Final radius (rho end for bobyqa, default: 0.001)",false,0.001,"optimizer final radius",cmd);
    TCLAP::ValueArg<double> searchStepArg("","st","Search step for exhaustive search (default: 2)",false,2,"exhaustive optimizer search step",cmd);

    TCLAP::ValueArg<double> translateUpperBoundArg("","tub","Upper bound on translation for bobyqa (in voxels, default: 10)",false,10,"Bobyqa translate upper bound",cmd);
    TCLAP::ValueArg<double> angleUpperBoundArg("","aub","Upper bound on angles for bobyqa (in degrees, default: 180)",false,180,"Bobyqa angle upper bound",cmd);
    TCLAP::ValueArg<double> scaleUpperBoundArg("","scu","Upper bound on scale for bobyqa (default: 3)",false,3,"Bobyqa scale upper bound",cmd);

    TCLAP::ValueArg<unsigned int> symmetryArg("","sym-reg","Registration symmetry type (0: asymmetric, 1: symmetric, 2: kissing, default: 0)",false,0,"symmetry type",cmd);

    TCLAP::ValueArg<unsigned int> agregatorArg("a","agregator","Transformation agregator type (0: M-Estimation, 1: least squares, 2: least trimmed squares, default: 0)",false,0,"agregator type",cmd);
    TCLAP::ValueArg<double> agregThresholdArg("","at","Agregator threshold value (for M-estimation or LTS)",false,0.5,"agregator threshold value",cmd);
    TCLAP::ValueArg<double> seStoppingThresholdArg("","lst","LTS Stopping Threshold (default: 0.01)",false,0.01,"LTS stopping threshold",cmd);

    TCLAP::ValueArg<unsigned int> numPyramidLevelsArg("p","pyr","Number of pyramid levels (default: 3)",false,3,"number of pyramid levels",cmd);
    TCLAP::ValueArg<unsigned int> lastPyramidLevelArg("l","last-level","Index of the last pyramid level explored (default: 0)",false,0,"last pyramid level",cmd);
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

    // Setting matcher arguments
    matcher->SetBlockSize( blockSizeArg.getValue() );
    matcher->SetBlockSpacing( blockSpacingArg.getValue() );
    matcher->SetStDevThreshold( stdevThresholdArg.getValue() );
    matcher->SetTransform( (PyramidBMType::Transform) blockTransfoArg.getValue() );
    matcher->SetAffineDirection(directionArg.getValue());
    matcher->SetMetric( (PyramidBMType::Metric) blockMetricArg.getValue() );
    matcher->SetOptimizer( (PyramidBMType::Optimizer) optimizerArg.getValue() );
    matcher->SetMaximumIterations( maxIterationsArg.getValue() );
    matcher->SetMinimalTransformError( minErrorArg.getValue() );
    matcher->SetFinalRadius(finalRadiusArg.getValue());
    matcher->SetOptimizerMaximumIterations( optimizerMaxIterationsArg.getValue() );
    matcher->SetSearchRadius( searchRadiusArg.getValue() );
    matcher->SetSearchAngleRadius( searchAngleRadiusArg.getValue() );
    matcher->SetSearchScaleRadius( searchScaleRadiusArg.getValue() );
    matcher->SetStepSize( searchStepArg.getValue() );
    matcher->SetTranslateUpperBound( translateUpperBoundArg.getValue() );
    matcher->SetAngleUpperBound( angleUpperBoundArg.getValue() );
    matcher->SetScaleUpperBound( scaleUpperBoundArg.getValue() );
    matcher->SetSymmetryType( (PyramidBMType::SymmetryType) symmetryArg.getValue() );
    matcher->SetAgregator( (PyramidBMType::Agregator) agregatorArg.getValue() );
    matcher->SetOutputTransformType( (PyramidBMType::OutputTransform) outTrTypeArg.getValue() );
    matcher->SetAgregThreshold( agregThresholdArg.getValue() );
    matcher->SetSeStoppingThreshold( seStoppingThresholdArg.getValue() );
    matcher->SetNumberOfPyramidLevels( numPyramidLevelsArg.getValue() );
    matcher->SetLastPyramidLevel( lastPyramidLevelArg.getValue() );

    if (blockMaskArg.getValue() != "")
        matcher->SetBlockGenerationMask(anima::readImage<PyramidBMType::MaskImageType>(blockMaskArg.getValue()));

    if (numThreadsArg.getValue() != 0)
        matcher->SetNumberOfThreads( numThreadsArg.getValue() );

    matcher->SetPercentageKept( percentageKeptArg.getValue() );
    matcher->SetInitializeOnCenterOfGravity( !initIdentityArg.isSet() );

    matcher->SetResultFile( outArg.getValue() );
    matcher->SetOutputTransformFile( outputTransformArg.getValue() );

    typedef itk::ImageFileReader<InputImageType> ReaderType;

    ReaderType::Pointer tmpRead = ReaderType::New();
    tmpRead->SetFileName(fixedArg.getValue());
    tmpRead->Update();

    matcher->SetReferenceImage(tmpRead->GetOutput());

    tmpRead = ReaderType::New();
    tmpRead->SetFileName(movingArg.getValue());
    tmpRead->Update();

    matcher->SetFloatingImage(tmpRead->GetOutput());

    if (initialTransformArg.getValue() != "")
        matcher->SetInitialTransform(initialTransformArg.getValue());

    AffineTransformPointer tmpTrsf = AffineTransformType::New();
    tmpTrsf->SetIdentity();

    matcher->SetOutputTransform(tmpTrsf.GetPointer());

    // Process
    itk::TimeProbe timer;
    timer.Start();

    try
    {
        matcher->Update();
        matcher->WriteOutputs();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }

    timer.Stop();

    std::cout << "Elapsed Time: " << timer.GetTotal()  << timer.GetUnit() << std::endl;

    return EXIT_SUCCESS;
}
