#include <animaPyramidalDenseMCMSVFMatchingBridge.h>
#include <animaMCMFileReader.h>
#include <animaGradientFileReader.h>

#include <tclap/CmdLine.h>
#include <itkTimeProbe.h>
#include <animaMCMConstants.h>

int main(int ac, const char** av)
{
    const unsigned int Dimension = 3;
    typedef anima::PyramidalDenseMCMSVFMatchingBridge <Dimension> PyramidBMType;
    typedef PyramidBMType::InputImageType InputImageType;
    typedef InputImageType::InternalPixelType InputInternalPixelType;

    // Parsing arguments
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    // Setting up parameters
    TCLAP::ValueArg<std::string> fixedArg("r","refimage","Fixed image",true,"","fixed image",cmd);
    TCLAP::ValueArg<std::string> movingArg("m","movingimage","Moving image",true,"","moving image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputimage","Output (registered) image",true,"","output image",cmd);
    TCLAP::ValueArg<std::string> outputTransformArg("O","outtransform","Output transformation",false,"","output transform",cmd);
    TCLAP::SwitchArg ppdImageArg("P","ppd", "Re-orientation strategy set to preservation of principal direction (default: finite strain)", cmd, false);

    TCLAP::ValueArg<unsigned int> blockSizeArg("","bs","Block size (default: 5)",false,5,"block size",cmd);
    TCLAP::ValueArg<unsigned int> blockSpacingArg("","sp","Block spacing (default: 2)",false,2,"block spacing",cmd);
    TCLAP::ValueArg<double> stdevThresholdArg("s","stdev","Threshold block standard deviation (default: 0.1)",false,0.1,"block minimal standard deviation",cmd);
    TCLAP::ValueArg<double> percentageKeptArg("k","per-kept","Percentage of blocks with the highest variance kept (default: 0.8)",false,0.8,"percentage of blocks kept",cmd);

    TCLAP::ValueArg<unsigned int> blockTransfoArg("t","in-transform","Transformation computed between blocks (0: translation, 1: rigid, 2: affine, default: 0)",false,0,"transformation between blocks",cmd);
    TCLAP::ValueArg<unsigned int> blockMetricArg("","metric","Similarity metric between blocks (0: basic mean squares, 1: one to on basic mean squares, 2: MCM mean squares, 3: MT pairing correlation, 4: MCM correlation, default: 1)",false,1,"similarity metric",cmd);
    TCLAP::ValueArg<unsigned int> blockOrientationArg("","bor","Re-orientation strategy when matching blocks (0: none, 1: finite strain, 2: PPD, default: 2)",false,2,"block re-orientation",cmd);

    TCLAP::ValueArg<double> smallDeltaArg("", "small-delta", "Diffusion small delta (in seconds)", false, anima::DiffusionSmallDelta, "small delta", cmd);
    TCLAP::ValueArg<double> bigDeltaArg("", "big-delta", "Diffusion big delta (in seconds)", false, anima::DiffusionBigDelta, "big delta", cmd);
    TCLAP::ValueArg<std::string> bvalArg("b","bvals","B-values",false,"","b-values",cmd);
    TCLAP::ValueArg<std::string> bvecArg("v","bvec","Gradient direction",false,"","gradient directions",cmd);

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
    TCLAP::ValueArg<double> neighborhoodApproximationArg("","na","Half size of the neighborhood approximation (multiplied by extrapolation sigma, default: 2.5)",false,2.5,"half size of neighborhood approximation",cmd);
    TCLAP::ValueArg<unsigned int> bchOrderArg("B","bch-order","BCH composition order (default: 1)",false,1,"BCH order",cmd);
    TCLAP::ValueArg<unsigned int> expOrderArg("e","exp-order","Order of field exponentiation approximation (in between 0 and 1, default: 0)",false,0,"exponentiation order",cmd);

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
        return EXIT_FAILURE;
    }

    PyramidBMType::Pointer matcher = PyramidBMType::New();

    anima::MCMFileReader <double,Dimension> refReader;
    refReader.SetFileName(fixedArg.getValue());
    refReader.Update();

    matcher->SetReferenceImage(refReader.GetModelVectorImage());

    anima::MCMFileReader <double,Dimension> floReader;
    floReader.SetFileName(movingArg.getValue());
    floReader.Update();

    matcher->SetFloatingImage(floReader.GetModelVectorImage());

    if ((bvalArg.getValue() != "")&&(bvecArg.getValue() != ""))
    {
        typedef anima::GradientFileReader <vnl_vector_fixed <double,3>, double> GradientReaderType;
        GradientReaderType gradientReader;
        gradientReader.SetGradientFileName(bvecArg.getValue());
        gradientReader.SetBValueBaseString(bvalArg.getValue());
        gradientReader.SetGradientIndependentNormalization(true);
        gradientReader.SetSmallDelta(smallDeltaArg.getValue());
        gradientReader.SetBigDelta(bigDeltaArg.getValue());

        gradientReader.Update();

        matcher->SetSmallDelta(smallDeltaArg.getValue());
        matcher->SetBigDelta(bigDeltaArg.getValue());
        matcher->SetGradientStrengths(gradientReader.GetGradientStrengths());
        matcher->SetGradientDirections(gradientReader.GetGradients());
    }

    // Setting matcher arguments
    matcher->SetBlockSize( blockSizeArg.getValue() );
    matcher->SetBlockSpacing( blockSpacingArg.getValue() );
    matcher->SetStDevThreshold( stdevThresholdArg.getValue() );
    matcher->SetTransform( (Transform) blockTransfoArg.getValue() );
    matcher->SetMetric( (Metric) blockMetricArg.getValue() );
    matcher->SetMetricOrientation( (MetricOrientationType) blockOrientationArg.getValue() );
    matcher->SetFiniteStrainImageReorientation(!ppdImageArg.isSet());
    matcher->SetOptimizer( (Optimizer) optimizerArg.getValue() );
    matcher->SetMaximumIterations( maxIterationsArg.getValue() );
    matcher->SetMinimalTransformError( minErrorArg.getValue() );
    matcher->SetOptimizerMaximumIterations( optimizerMaxIterationsArg.getValue() );
    matcher->SetStepSize( searchStepArg.getValue() );
    matcher->SetTranslateUpperBound( translateUpperBoundArg.getValue() );
    matcher->SetAngleUpperBound( angleUpperBoundArg.getValue() );
    matcher->SetScaleUpperBound( scaleUpperBoundArg.getValue() );
    matcher->SetSymmetryType( (SymmetryType) symmetryArg.getValue() );
    matcher->SetAgregator( (Agregator) agregatorArg.getValue() );
    matcher->SetExtrapolationSigma(extrapolationSigmaArg.getValue());
    matcher->SetElasticSigma(elasticSigmaArg.getValue());
    matcher->SetOutlierSigma(outlierSigmaArg.getValue());
    matcher->SetMEstimateConvergenceThreshold(mEstimateConvergenceThresholdArg.getValue());
    matcher->SetNeighborhoodApproximation(neighborhoodApproximationArg.getValue());
    matcher->SetBCHCompositionOrder(bchOrderArg.getValue());
    matcher->SetExponentiationOrder(expOrderArg.getValue());
    matcher->SetNumberOfPyramidLevels( numPyramidLevelsArg.getValue() );
    matcher->SetLastPyramidLevel( lastPyramidLevelArg.getValue() );

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
