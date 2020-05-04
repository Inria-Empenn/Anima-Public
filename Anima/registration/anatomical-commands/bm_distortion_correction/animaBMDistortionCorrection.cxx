#include <tclap/CmdLine.h>

#include <animaPyramidalDistortionCorrectionBlockMatchingBridge.h>
#include <animaReadWriteFunctions.h>

#include <itkTimeProbe.h>

int main(int argc, const char** argv)
{
    typedef anima::PyramidalDistortionCorrectionBlockMatchingBridge <3> PyramidBMType;
    typedef PyramidBMType::InputImageType InputImageType;
    typedef PyramidBMType::VectorFieldType VectorFieldType;

    // Parsing arguments
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    
    // Setting up parameters
    TCLAP::ValueArg<std::string> backwardArg("b","backimage","Backward image (eg PA image)",true,"","backward image",cmd);
    TCLAP::ValueArg<std::string> forwardArg("f","forwardimage","Forward image (eg AP image)",true,"","forward image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputimage","Output (corrected) image",true,"","output image",cmd);
    TCLAP::ValueArg<std::string> initialTransformArg("i","initransform","Initial transformation",false,"","initial transform",cmd);
    TCLAP::ValueArg<std::string> outputTransformArg("O","outtransform","Output transformation",false,"","output transform",cmd);
    
    TCLAP::ValueArg<unsigned int> blockSizeArg("","bs","Block size (default: 3)",false,3,"block size",cmd);
    TCLAP::ValueArg<unsigned int> blockSpacingArg("","sp","Block spacing (default: 2)",false,2,"block spacing",cmd);
    TCLAP::ValueArg<double> stdevThresholdArg("s","stdev","Threshold block standard deviation (default: 15)",false,15,"block minimal standard deviation",cmd);
    TCLAP::ValueArg<double> percentageKeptArg("k","per-kept","Percentage of blocks with the highest variance kept (default: 0.8)",false,0.8,"percentage of blocks kept",cmd);
    
    TCLAP::ValueArg<unsigned int> maxIterationsArg("","mi","Maximum block match iterations (default: 10)",false,10,"maximum iterations",cmd);
    
    TCLAP::ValueArg<unsigned int> directionArg("d","dir","Gradient phase direction in image (default: 1 = Y axis)",false,1,"gradient phase direction",cmd);
    TCLAP::ValueArg<unsigned int> optimizerMaxIterationsArg("","oi","Maximum iterations for local optimizer (default: 100)",false,100,"maximum local optimizer iterations",cmd);
    
    TCLAP::ValueArg<double> searchRadiusArg("","sr","Search radius in pixels (rho start for bobyqa, default: 2)",false,2,"optimizer search radius",cmd);
    TCLAP::ValueArg<double> searchScaleRadiusArg("","scr","Search scale radius (rho start for bobyqa, default: 0.1)",false,0.1,"optimizer search scale radius",cmd);
    TCLAP::ValueArg<double> searchSkewRadiusArg("","skr","Search skew radius (rho start for bobyqa, default: 0.1)",false,0.1,"optimizer search skew radius",cmd);
    TCLAP::ValueArg<double> finalRadiusArg("","fr","Final radius (rho end for bobyqa, default: 0.001)",false,0.001,"optimizer final radius",cmd);
    
    TCLAP::ValueArg<double> translateUpperBoundArg("","tub","Upper bound on translation for bobyqa (in voxels, default: 10)",false,10,"Bobyqa translate upper bound",cmd);
    TCLAP::ValueArg<double> scaleUpperBoundArg("","scu","Upper bound on scale for bobyqa (default: 5)",false,5,"Bobyqa scale upper bound",cmd);
    TCLAP::ValueArg<double> skewUpperBoundArg("","sku","Upper bound on skew for bobyqa (in degrees, default: 45)",false,45,"Bobyqa skew upper bound",cmd);

    TCLAP::ValueArg<unsigned int> blockTransfoArg("t","in-transform","Transformation computed between blocks (0: direction, 1: direction+scale, 2: direction+scale+skew, default: 2)",false,2,"transformation between blocks",cmd);
    TCLAP::ValueArg<unsigned int> agregatorArg("","agregator","Transformation agregator type (0: Baloo, 1: M-smoother, default: 0)",false,0,"agregator type",cmd);
    TCLAP::ValueArg<unsigned int> blockMetricArg("","metric","Similarity metric between blocks (0: correlation coefficient, 1: squared correlation coefficient, 2: mean squares, default: 1)",false,1,"similarity metric",cmd);

    TCLAP::SwitchArg weightedAgregationArg("w","no-weighted-agregation", "If set, weighted agregation is deactivated", cmd, false);
    TCLAP::ValueArg<double> extrapolationSigmaArg("","fs","Sigma for extrapolation of local pairings (default: 3)",false,3,"extrapolation sigma",cmd);
    TCLAP::ValueArg<double> elasticSigmaArg("","es","Sigma for elastic regularization (default: 2)",false,2,"elastic regularization sigma",cmd);
    TCLAP::ValueArg<double> outlierSigmaArg("","os","Sigma for outlier rejection among local pairings (default: 3)",false,3,"outlier rejection sigma",cmd);
    TCLAP::ValueArg<double> mEstimateConvergenceThresholdArg("","met","Threshold to consider m-estimator converged (default: 0.01)",false,0.01,"m-estimation convergence threshold",cmd);
    TCLAP::ValueArg<double> neighborhoodApproximationArg("","na","Half size of the neighborhood approximation (multiplied by extrapolation sigma, default: 2.5)",false,2.5,"half size of neighborhood approximation",cmd);
    
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

    PyramidBMType *matcher = new PyramidBMType;

    matcher->SetBackwardImage(anima::readImage<InputImageType>(backwardArg.getValue()).GetPointer());
    matcher->SetForwardImage(anima::readImage<InputImageType>(forwardArg.getValue()).GetPointer());
    
    if (initialTransformArg.getValue() != "")
        matcher->SetInitialTransformField(anima::readImage<VectorFieldType>(initialTransformArg.getValue()));
    
    // Setting matcher arguments
    matcher->SetBlockSize( blockSizeArg.getValue() );
    matcher->SetBlockSpacing( blockSpacingArg.getValue() );
    matcher->SetStDevThreshold( stdevThresholdArg.getValue() );
    matcher->SetMaximumIterations( maxIterationsArg.getValue() );
    matcher->SetOptimizerMaximumIterations( optimizerMaxIterationsArg.getValue() );
    matcher->SetSearchRadius( searchRadiusArg.getValue() );
    matcher->SetSearchScaleRadius( searchScaleRadiusArg.getValue() );
    matcher->SetSearchSkewRadius( searchSkewRadiusArg.getValue() );
    matcher->SetFinalRadius(finalRadiusArg.getValue());
    matcher->SetTranlateUpperBound( translateUpperBoundArg.getValue() );
    matcher->SetScaleUpperBound( std::log(scaleUpperBoundArg.getValue()) );
    matcher->SetSkewUpperBound( std::tan(skewUpperBoundArg.getValue() * M_PI / 180.0) );

    matcher->SetTransformKind( (TransformKind) blockTransfoArg.getValue() );
    matcher->SetAgregator( (Agregator) agregatorArg.getValue() );
    matcher->SetMetric( (Metric) blockMetricArg.getValue() );
    
    matcher->SetExtrapolationSigma(extrapolationSigmaArg.getValue());
    matcher->SetElasticSigma(elasticSigmaArg.getValue());
    matcher->SetOutlierSigma(outlierSigmaArg.getValue());
    matcher->SetMEstimateConvergenceThreshold(mEstimateConvergenceThresholdArg.getValue());
    matcher->SetNeighborhoodApproximation(neighborhoodApproximationArg.getValue());
    matcher->SetWeightedAgregation( weightedAgregationArg.isSet() );
    matcher->SetNumberOfPyramidLevels( numPyramidLevelsArg.getValue() );
    matcher->SetLastPyramidLevel( lastPyramidLevelArg.getValue() );
    matcher->SetTransformDirection(directionArg.getValue());
    matcher->SetExponentiationOrder(expOrderArg.getValue());

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
    catch(itk::ExceptionObject &e)
    {
        std::cerr << "Error " << e << std::endl;
        return EXIT_FAILURE;
    }

    timer.Stop();

    std::cout << "Elapsed Time: " << timer.GetTotal()  << timer.GetUnit() << std::endl;

    return EXIT_SUCCESS;
}
