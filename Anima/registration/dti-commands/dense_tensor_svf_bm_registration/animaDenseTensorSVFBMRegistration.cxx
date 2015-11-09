#include <animaPyramidalDenseTensorSVFMatchingBridge.h>

#include <tclap/CmdLine.h>

#include <itkTimeProbe.h>

int main(int ac, const char** av)
{
    const    unsigned int    Dimension = 3;
    typedef anima::PyramidalDenseTensorSVFMatchingBridge <Dimension> PyramidBMType;
    typedef PyramidBMType::InputImageType InputImageType;
    typedef InputImageType::InternalPixelType InputInternalPixelType;
    typedef itk::ImageFileReader<InputImageType> ReaderType;

    // Parsing arguments
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',ANIMA_VERSION);

    // Setting up parameters
    TCLAP::ValueArg<std::string> fixedArg("r","refimage","Fixed image",true,"","fixed image",cmd);
    TCLAP::ValueArg<std::string> movingArg("m","movingimage","Moving image",true,"","moving image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputimage","Output (registered) image",true,"","output image",cmd);
    TCLAP::ValueArg<std::string> outputTransformArg("O","outtransform","Output transformation",false,"","output transform",cmd);

    TCLAP::ValueArg<unsigned int> blockSizeArg("","bs","Block size (default: 5)",false,5,"block size",cmd);
    TCLAP::ValueArg<unsigned int> blockSpacingArg("","sp","Block spacing (default: 2)",false,2,"block spacing",cmd);
    TCLAP::ValueArg<float> stdevThresholdArg("s","stdev","Threshold block standard deviation (default: 0)",false,0,"block minimal standard deviation",cmd);
    TCLAP::ValueArg<double> percentageKeptArg("k","per-kept","Percentage of blocks with the highest variance kept (default: 0.8)",false,0.8,"percentage of blocks kept",cmd);

    TCLAP::ValueArg<unsigned int> blockTransfoArg("t","in-transform","Transformation computed between blocks (0: translation, 1: rigid, 2: affine, default: 0)",false,0,"transformation between blocks",cmd);
    TCLAP::ValueArg<unsigned int> blockMetricArg("","metric","Similarity metric between blocks (0: mean squares, 1: tensor correlation, 2: tensor generalized correlation, 3: tensor oriented generalized correlation, default: 3)",false,3,"similarity metric",cmd);
    TCLAP::ValueArg<unsigned int> optimizerArg("","opt","Optimizer for optimal block search (0: Exhaustive, 1: Bobyqa, default: 1)",false,1,"optimizer",cmd);

    TCLAP::ValueArg<unsigned int> maxIterationsArg("","mi","Maximum block match iterations (default: 10)",false,10,"maximum iterations",cmd);
    TCLAP::ValueArg<float> minErrorArg("","me","Minimal distance between consecutive estimated transforms (default: 0.01)",false,0.01,"minimal distance between transforms",cmd);

    TCLAP::ValueArg<unsigned int> optimizerMaxIterationsArg("","oi","Maximum iterations for local optimizer (default: 100)",false,100,"maximum local optimizer iterations",cmd);

    TCLAP::ValueArg<double> searchRadiusArg("","sr","Search radius in pixels (exhaustive search window, rho start for bobyqa, default: 2)",false,2,"optimizer search radius",cmd);
    TCLAP::ValueArg<double> searchAngleRadiusArg("","sar","Search angle radius in degrees (rho start for bobyqa, default: 5)",false,5,"optimizer search angle radius",cmd);
    TCLAP::ValueArg<double> searchScaleRadiusArg("","scr","Search scale radius (rho start for bobyqa, default: 0.1)",false,0.1,"optimizer search scale radius",cmd);
    TCLAP::ValueArg<double> searchSkewRadiusArg("","skr","Search skew radius in degrees (rho start for bobyqa, default: 5)",false,5,"optimizer search skew radius",cmd);
    TCLAP::ValueArg<double> finalRadiusArg("","fr","Final radius (rho end for bobyqa, default: 0.001)",false,0.001,"optimizer final radius",cmd);
    TCLAP::ValueArg<double> searchStepArg("","st","Search step for exhaustive search (default: 1)",false,1,"exhaustive optimizer search step",cmd);

    TCLAP::ValueArg<double> translateUpperBoundArg("","tub","Upper bound on translation for bobyqa (in voxels, default: 10)",false,10,"Bobyqa translate upper bound",cmd);
    TCLAP::ValueArg<double> angleUpperBoundArg("","aub","Upper bound on angles for bobyqa (in degrees, default: 180)",false,180,"Bobyqa angle upper bound",cmd);
    TCLAP::ValueArg<double> scaleUpperBoundArg("","scu","Upper bound on scale for bobyqa (default: 3)",false,3,"Bobyqa scale upper bound",cmd);
    TCLAP::ValueArg<double> skewUpperBoundArg("","sku","Upper bound on skew for bobyqa (in degrees, default: 45)",false,45,"Bobyqa skew upper bound",cmd);
    
    TCLAP::ValueArg<unsigned int> symmetryArg("","sym-reg","Registration symmetry type (0: asymmetric, 1: symmetric, 2: kissing, default: 0)",false,0,"symmetry type",cmd);
    TCLAP::ValueArg<unsigned int> agregatorArg("a","agregator","Transformation agregator type (0: Baloo, 1: M-smoother, default: 0)",false,0,"agregator type",cmd);
    TCLAP::ValueArg<double> extrapolationSigmaArg("","fs","Sigma for extrapolation of local pairings (default: 3)",false,3,"extrapolation sigma",cmd);
    TCLAP::ValueArg<double> elasticSigmaArg("","es","Sigma for elastic regularization (default: 3)",false,3,"elastic regularization sigma",cmd);
    TCLAP::ValueArg<double> outlierSigmaArg("","os","Sigma for outlier rejection among local pairings (default: 3)",false,3,"outlier rejection sigma",cmd);
    TCLAP::ValueArg<double> mEstimateConvergenceThresholdArg("","met","Threshold to consider m-estimator converged (default: 0.01)",false,0.01,"m-estimation convergence threshold",cmd);
    TCLAP::ValueArg<double> neighborhoodApproximationArg("","na","Half size of the neighborhood approximation (multiplied by extrapolation sigma, default: 2.5)",false,2.5,"half size of neighborhood approximation",cmd);

    TCLAP::SwitchArg useTransformDamArg("D","use-dam", "Activate transformation dam to force identity far away from any blocks", cmd, false);
    TCLAP::ValueArg<double> damDistanceArg("","dd","Distance of the deformation dam from the outer blocks (will be multiplied by extrapolation sigma, default: 2.5)",false,2.5,"identity dam distance",cmd);

    TCLAP::ValueArg<unsigned int> numPyramidLevelsArg("p","pyr","Number of pyramid levels (default: 3)",false,3,"number of pyramid levels",cmd);
    TCLAP::ValueArg<unsigned int> lastPyramidLevelArg("l","last-level","Index of the last pyramid level explored (default: 0)",false,0,"last pyramid level",cmd);
    TCLAP::ValueArg<unsigned int> numThreadsArg("T","threads","Number of execution threads (default: 0 = all cores)",false,0,"number of threads",cmd);

    try
    {
        cmd.parse(ac,av);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        exit(-1);
    }

    PyramidBMType::Pointer matcher = PyramidBMType::New();

    ReaderType::Pointer tmpRead = ReaderType::New();
    tmpRead->SetFileName(fixedArg.getValue());
    tmpRead->Update();

    typedef anima::LogTensorImageFilter <InputInternalPixelType,Dimension> LogTensorFilterType;
    LogTensorFilterType::Pointer tensorRefLogger = LogTensorFilterType::New();

    tensorRefLogger->SetInput(tmpRead->GetOutput());
    tensorRefLogger->SetScaleNonDiagonal(true);

    if (numThreadsArg.getValue() != 0)
        tensorRefLogger->SetNumberOfThreads(numThreadsArg.getValue());

    tensorRefLogger->Update();

    matcher->SetReferenceImage(tensorRefLogger->GetOutput());
    matcher->GetReferenceImage()->DisconnectPipeline();

    tmpRead = ReaderType::New();
    tmpRead->SetFileName(movingArg.getValue());
    tmpRead->Update();

    LogTensorFilterType::Pointer tensorFloLogger = LogTensorFilterType::New();

    tensorFloLogger->SetInput(tmpRead->GetOutput());
    tensorFloLogger->SetScaleNonDiagonal(true);

    if (numThreadsArg.getValue() != 0)
        tensorFloLogger->SetNumberOfThreads(numThreadsArg.getValue());

    tensorFloLogger->Update();

    matcher->SetFloatingImage(tensorFloLogger->GetOutput());
    matcher->GetFloatingImage()->DisconnectPipeline();

    // Setting matcher arguments
    matcher->SetBlockSize( blockSizeArg.getValue() );
    matcher->SetBlockSpacing( blockSpacingArg.getValue() );
    matcher->SetStDevThreshold( stdevThresholdArg.getValue() );
    matcher->SetTransform( (Transform) blockTransfoArg.getValue() );
    matcher->SetMetric( (Metric) blockMetricArg.getValue() );
    matcher->SetOptimizer( (Optimizer) optimizerArg.getValue() );
    matcher->SetMaximumIterations( maxIterationsArg.getValue() );
    matcher->SetMinimalTransformError( minErrorArg.getValue() );
    matcher->SetFinalRadius(finalRadiusArg.getValue());
    matcher->SetOptimizerMaximumIterations( optimizerMaxIterationsArg.getValue() );
    matcher->SetSearchRadius( searchRadiusArg.getValue() );
    matcher->SetSearchAngleRadius( searchAngleRadiusArg.getValue() );
    matcher->SetSearchSkewRadius( searchSkewRadiusArg.getValue() );
    matcher->SetSearchScaleRadius( searchScaleRadiusArg.getValue() );
    matcher->SetStepSize( searchStepArg.getValue() );
    matcher->SetTranslateUpperBound( translateUpperBoundArg.getValue() );
    matcher->SetAngleUpperBound( angleUpperBoundArg.getValue() );
    matcher->SetSkewUpperBound( skewUpperBoundArg.getValue() );
    matcher->SetScaleUpperBound( scaleUpperBoundArg.getValue() );
    matcher->SetSymmetryType( (SymmetryType) symmetryArg.getValue() );
    matcher->SetAgregator( (Agregator) agregatorArg.getValue() );
    matcher->SetExtrapolationSigma(extrapolationSigmaArg.getValue());
    matcher->SetElasticSigma(elasticSigmaArg.getValue());
    matcher->SetOutlierSigma(outlierSigmaArg.getValue());
    matcher->SetMEstimateConvergenceThreshold(mEstimateConvergenceThresholdArg.getValue());
    matcher->SetNeighborhoodApproximation(neighborhoodApproximationArg.getValue());
    matcher->SetUseTransformationDam(useTransformDamArg.isSet());
    matcher->SetDamDistance(damDistanceArg.getValue());
    matcher->SetNumberOfPyramidLevels( numPyramidLevelsArg.getValue() );
    matcher->SetLastPyramidLevel( lastPyramidLevelArg.getValue() );

    if (numThreadsArg.getValue() != 0)
        matcher->SetNumberOfThreads( numThreadsArg.getValue() );

    matcher->SetPercentageKept( percentageKeptArg.getValue() );

    matcher->SetResultFile( outArg.getValue() );
    matcher->SetOutputTransformFile( outputTransformArg.getValue() );

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
        return 1;
    }

    timer.Stop();

    std::cout << "Elapsed Time: " << timer.GetTotal()  << timer.GetUnit() << std::endl;

    return EXIT_SUCCESS;
}
