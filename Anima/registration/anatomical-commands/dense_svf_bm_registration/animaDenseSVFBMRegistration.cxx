#include <animaPyramidalDenseSVFMatchingBridge.h>

#include <tclap/CmdLine.h>

#include <itkTimeProbe.h>

int main(int argc, const char** argv)
{
    typedef anima::PyramidalDenseSVFMatchingBridge <3> PyramidBMType;
    typedef PyramidBMType::InputImageType InputImageType;
    typedef itk::ImageFileReader<InputImageType> ReaderType;

    // Parsing arguments
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    // Setting up parameters
    TCLAP::ValueArg<std::string> fixedArg("r","refimage","Fixed image",true,"","fixed image",cmd);
    TCLAP::ValueArg<std::string> movingArg("m","movingimage","Moving image",true,"","moving image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputimage","Output (registered) image",true,"","output image",cmd);
    TCLAP::ValueArg<std::string> outputTransformArg("O","outtransform","Output transformation",false,"","output transform",cmd);
    TCLAP::ValueArg<std::string> blockMaskArg("M","mask-im","Mask image for block generation",false,"","block mask image",cmd);

    TCLAP::ValueArg<unsigned int> blockSizeArg("","bs","Block size (default: 5)",false,5,"block size",cmd);
    TCLAP::ValueArg<unsigned int> blockSpacingArg("","sp","Block spacing (default: 2)",false,2,"block spacing",cmd);
    TCLAP::ValueArg<double> stdevThresholdArg("s","stdev","Threshold block standard deviation (default: 5)",false,5,"block minimal standard deviation",cmd);
    TCLAP::ValueArg<double> percentageKeptArg("k","per-kept","Percentage of blocks with the highest variance kept (default: 0.8)",false,0.8,"percentage of blocks kept",cmd);

    TCLAP::ValueArg<unsigned int> blockTransfoArg("t","in-transform","Transformation computed between blocks (0: translation, 1: rigid, 2: affine, 3: directional affine, default: 0)",false,0,"transformation between blocks",cmd);
    TCLAP::ValueArg<unsigned int> directionArg("d","dir","Affine direction for directional transform output (default: 1 = Y axis)",false,1,"direction of directional affine",cmd);
    TCLAP::ValueArg<unsigned int> blockMetricArg("","metric","Similarity metric between blocks (0: mean squares, 1: correlation coefficient, 2: squared correlation coefficient, default: 2)",false,2,"similarity metric",cmd);
    TCLAP::ValueArg<unsigned int> optimizerArg("","opt","Optimizer for optimal block search (0: Exhaustive, 1: Bobyqa, default: 1)",false,1,"optimizer",cmd);

    TCLAP::ValueArg<unsigned int> maxIterationsArg("","mi","Maximum block match iterations (default: 10)",false,10,"maximum iterations",cmd);
    TCLAP::ValueArg<double> minErrorArg("","me","Minimal distance between consecutive estimated transforms (default: 0.01)",false,0.01,"minimal distance between transforms",cmd);

    TCLAP::ValueArg<unsigned int> optimizerMaxIterationsArg("","oi","Maximum iterations for local optimizer (default: 100)",false,100,"maximum local optimizer iterations",cmd);

    TCLAP::ValueArg<double> searchStepArg("","st","Search step for exhaustive search (default: 1)",false,1,"exhaustive optimizer search step",cmd);

    TCLAP::ValueArg<double> translateUpperBoundArg("","tub","Upper bound on translation for bobyqa (in voxels, default: 3)",false,3,"Bobyqa translate upper bound",cmd);
    TCLAP::ValueArg<double> angleUpperBoundArg("","aub","Upper bound on angles for bobyqa (in degrees, default: 180)",false,180,"Bobyqa angle upper bound",cmd);
    TCLAP::ValueArg<double> scaleUpperBoundArg("","scu","Upper bound on scale for bobyqa (default: 3)",false,3,"Bobyqa scale upper bound",cmd);

    TCLAP::ValueArg<unsigned int> symmetryArg("","sym-reg","Registration symmetry type (0: asymmetric, 1: symmetric, 2: kissing, default: 0)",false,0,"symmetry type",cmd);
    TCLAP::ValueArg<double> kissingLocationArg("K","kissing-point","Kissing point location along the transformation path (default: half-way = 0.5)",false,0.5,"kissing point location",cmd);

    TCLAP::ValueArg<unsigned int> agregatorArg("a","agregator","Transformation agregator type (0: Baloo, 1: M-smoother, default: 0)",false,0,"agregator type",cmd);
    TCLAP::ValueArg<double> extrapolationSigmaArg("","fs","Sigma for extrapolation of local pairings (default: 3)",false,3,"extrapolation sigma",cmd);
    TCLAP::ValueArg<double> elasticSigmaArg("","es","Sigma for elastic regularization (default: 3)",false,3,"elastic regularization sigma",cmd);
    TCLAP::ValueArg<double> outlierSigmaArg("","os","Sigma for outlier rejection among local pairings (default: 3)",false,3,"outlier rejection sigma",cmd);
    TCLAP::ValueArg<double> mEstimateConvergenceThresholdArg("","met","Threshold to consider m-estimator converged (default: 0.01)",false,0.01,"m-estimation convergence threshold",cmd);
    TCLAP::ValueArg<unsigned int> bchOrderArg("b","bch-order","BCH composition order (default: 1)",false,1,"BCH order",cmd);
    TCLAP::ValueArg<unsigned int> expOrderArg("e","exp-order","Order of field exponentiation approximation (in between 0 and 1, default: 0)",false,0,"exponentiation order",cmd);

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

    PyramidBMType::Pointer matcher = PyramidBMType::New();

    ReaderType::Pointer tmpRead = ReaderType::New();
    tmpRead->SetFileName(fixedArg.getValue());
    tmpRead->Update();

    matcher->SetReferenceImage(tmpRead->GetOutput());

    tmpRead = ReaderType::New();
    tmpRead->SetFileName(movingArg.getValue());
    tmpRead->Update();

    matcher->SetFloatingImage(tmpRead->GetOutput());

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
    matcher->SetOptimizerMaximumIterations( optimizerMaxIterationsArg.getValue() );
    matcher->SetStepSize( searchStepArg.getValue() );
    matcher->SetTranslateUpperBound( translateUpperBoundArg.getValue() );
    matcher->SetAngleUpperBound( angleUpperBoundArg.getValue() );
    matcher->SetScaleUpperBound( scaleUpperBoundArg.getValue() );
    matcher->SetSymmetryType( (PyramidBMType::SymmetryType) symmetryArg.getValue() );
    matcher->SetAgregator( (PyramidBMType::Agregator) agregatorArg.getValue() );
    matcher->SetExtrapolationSigma(extrapolationSigmaArg.getValue());
    matcher->SetElasticSigma(elasticSigmaArg.getValue());
    matcher->SetOutlierSigma(outlierSigmaArg.getValue());
    matcher->SetMEstimateConvergenceThreshold(mEstimateConvergenceThresholdArg.getValue());
    matcher->SetBCHCompositionOrder(bchOrderArg.getValue());
    matcher->SetExponentiationOrder(expOrderArg.getValue());
    matcher->SetNumberOfPyramidLevels( numPyramidLevelsArg.getValue() );
    matcher->SetLastPyramidLevel( lastPyramidLevelArg.getValue() );
    matcher->SetRegistrationPointLocation(kissingLocationArg.getValue());

    if (blockMaskArg.getValue() != "")
        matcher->SetBlockGenerationMask(anima::readImage<PyramidBMType::MaskImageType>(blockMaskArg.getValue()));

    if (numThreadsArg.getValue() != 0)
        matcher->SetNumberOfWorkUnits( numThreadsArg.getValue() );

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
        return EXIT_FAILURE;
    }

    timer.Stop();

    std::cout << "Elapsed Time: " << timer.GetTotal()  << timer.GetUnit() << std::endl;

    return EXIT_SUCCESS;
}
