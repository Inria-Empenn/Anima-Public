#include <tclap/CmdLine.h>

#include <animaPyramidalBlockMatchingBridge.h>
#include <animaPyramidalDenseSVFMatchingBridge.h>
#include <animaReadWriteFunctions.h>
#include <itkExtractImageFilter.h>

#include <itkImageRegionIterator.h>
#include <itkCompositeTransform.h>
#include <itkStationaryVelocityFieldTransform.h>
#include <rpiDisplacementFieldTransform.h>
#include <animaVelocityUtils.h>
#include <animaResampleImageFilter.h>
#include <animaGradientFileReader.h>

int main(int argc, const char** argv)
{
    const unsigned int Dimension = 3;

    typedef itk::Image <double,Dimension+1> InputImageType;
    typedef itk::Image <double,Dimension> InputSubImageType;
    typedef itk::ImageRegionIterator <InputImageType> InputImageIteratorType;
    typedef itk::ImageRegionIterator <InputSubImageType> InputSubImageIteratorType;

    typedef anima::PyramidalBlockMatchingBridge <Dimension> PyramidBMType;
    typedef anima::BaseTransformAgregator <Dimension> AgregatorType;
    typedef itk::AffineTransform<AgregatorType::ScalarType,Dimension> AffineTransformType;
    typedef AffineTransformType::Pointer AffineTransformPointer;

    // Parsing arguments
    TCLAP::CmdLine  cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    // Setting up parameters
    TCLAP::ValueArg<std::string> inputArg("i","input","Input 4D image",true,"","input 4D image",cmd);
    TCLAP::ValueArg<std::string> inBVecArg("I","input-bvec","Input gradient vectors file",true,"","input gradients",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","Output (corrected) image",true,"","output image",cmd);
    TCLAP::ValueArg<std::string> outBVecArg("O","output-bvec","Output gradient vectors (bvec format)",true,"","output gradients",cmd);

    TCLAP::ValueArg<unsigned int> directionArg("d","dir","Affine direction for directional transform output (default: 1 = Y axis)",false,1,"direction of directional affine",cmd);
    TCLAP::ValueArg<unsigned int> b0Arg("b","b0","Index of the B0 reference image",false,0,"reference image index",cmd);

    TCLAP::ValueArg<unsigned int> blockSizeArg("","bs","Block size (default: 5)",false,5,"block size",cmd);
    TCLAP::ValueArg<unsigned int> blockSpacingArg("","sp","Block spacing (default: 5)",false,5,"block spacing",cmd);
    TCLAP::ValueArg<unsigned int> nlBlockSpacingArg("","nsp","Block spacing (default: 3)",false,3,"non linear matching block spacing",cmd);
    TCLAP::ValueArg<double> stdevThresholdArg("s","stdev","Threshold block standard deviation (default: 5)",false,5,"block minimal standard deviation",cmd);
    TCLAP::ValueArg<double> percentageKeptArg("k","per-kept","Percentage of blocks with the highest variance kept (default: 0.8)",false,0.8,"percentage of blocks kept",cmd);

    TCLAP::ValueArg<unsigned int> blockMetricArg("","metric","Similarity metric between blocks (0: mean squares, 1: correlation coefficient, 2: squared correlation coefficient, default: 2)",false,2,"similarity metric",cmd);
    TCLAP::ValueArg<unsigned int> optimizerArg("","opt","Optimizer for optimal block search (0: Exhaustive, 1: Bobyqa, default: 1)",false,1,"optimizer",cmd);

    TCLAP::ValueArg<unsigned int> maxIterationsArg("","mi","Maximum block match iterations (default: 10)",false,10,"maximum iterations",cmd);
    TCLAP::ValueArg<double> minErrorArg("","me","Minimal distance between consecutive estimated transforms (default: 0.01)",false,0.01,"minimal distance between transforms",cmd);

    TCLAP::ValueArg<unsigned int> optimizerMaxIterationsArg("","oi","Maximum iterations for local optimizer (default: 100)",false,100,"maximum local optimizer iterations",cmd);

    TCLAP::ValueArg<double> searchStepArg("","st","Search step for exhaustive search (default: 2)",false,2,"exhaustive optimizer search step",cmd);
    TCLAP::ValueArg<double> translateUpperBoundArg("","tub","Upper bound on translation for bobyqa (in voxels, default: 3)",false,3,"Bobyqa translate upper bound",cmd);

    TCLAP::ValueArg<unsigned int> symmetryArg("","sym-reg","Registration symmetry type (0: asymmetric, 1: symmetric, 2: kissing, default: 0)",false,0,"symmetry type",cmd);
    TCLAP::ValueArg<unsigned int> agregatorArg("a","agregator","Transformation agregator type (0: M-Estimation, 1: least squares, 2: least trimmed squares, default: 0)",false,0,"agregator type",cmd);
    TCLAP::ValueArg<double> agregThresholdArg("","at","Agregator threshold value (for M-estimation or LTS)",false,0.5,"agregator threshold value",cmd);
    TCLAP::ValueArg<double> extrapolationSigmaArg("","fs","Sigma for extrapolation of local pairings (default: 3)",false,3,"extrapolation sigma",cmd);
    TCLAP::ValueArg<double> elasticSigmaArg("","es","Sigma for elastic regularization (default: 3)",false,3,"elastic regularization sigma",cmd);
    TCLAP::ValueArg<double> outlierSigmaArg("","os","Sigma for outlier rejection among local pairings (default: 3)",false,3,"outlier rejection sigma",cmd);
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

    InputImageType::Pointer inputImage = anima::readImage <InputImageType> (inputArg.getValue());
    unsigned int numberOfImages = inputImage->GetLargestPossibleRegion().GetSize()[Dimension];
    typedef itk::ExtractImageFilter <InputImageType, InputSubImageType> ExtractFilterType;

    ExtractFilterType::Pointer referenceExtractFilter = ExtractFilterType::New();
    referenceExtractFilter->SetInput(inputImage);
    InputImageType::RegionType refExtractRegion = inputImage->GetLargestPossibleRegion();
    refExtractRegion.SetIndex(Dimension,b0Arg.getValue());
    refExtractRegion.SetSize(Dimension,0);
    referenceExtractFilter->SetExtractionRegion(refExtractRegion);
    referenceExtractFilter->SetDirectionCollapseToGuess();
    referenceExtractFilter->Update();

    typedef anima::GradientFileReader < vnl_vector_fixed <double,3>, double > GFReaderType;
    GFReaderType gfReader;
    gfReader.SetGradientFileName(inBVecArg.getValue());
    gfReader.SetGradientIndependentNormalization(false);
    gfReader.Update();

    GFReaderType::GradientVectorType directions = gfReader.GetGradients();

    for (unsigned int i = 0;i < numberOfImages;++i)
    {
        if (i == b0Arg.getValue())
            continue;

        std::cout << "\033[K\rProcessing image " << i+1 << " out of " << numberOfImages << std::flush;

        ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
        extractFilter->SetInput(inputImage);
        InputImageType::RegionType extractRegion = inputImage->GetLargestPossibleRegion();
        extractRegion.SetIndex(Dimension,i);
        extractRegion.SetSize(Dimension,0);
        extractFilter->SetExtractionRegion(extractRegion);
        extractFilter->SetDirectionCollapseToGuess();
        extractFilter->Update();

        // First perform rigid registration to correct for movement
        PyramidBMType::Pointer matcher = PyramidBMType::New();

        matcher->SetBlockSize(blockSizeArg.getValue());
        matcher->SetBlockSpacing(blockSpacingArg.getValue());
        matcher->SetStDevThreshold(stdevThresholdArg.getValue());
        matcher->SetMetric((PyramidBMType::Metric) blockMetricArg.getValue());
        matcher->SetOptimizer((PyramidBMType::Optimizer) optimizerArg.getValue());
        matcher->SetMaximumIterations(maxIterationsArg.getValue());
        matcher->SetMinimalTransformError(minErrorArg.getValue());
        matcher->SetOptimizerMaximumIterations(optimizerMaxIterationsArg.getValue());
        matcher->SetStepSize(searchStepArg.getValue());
        matcher->SetTranslateUpperBound(translateUpperBoundArg.getValue());
        matcher->SetSymmetryType((PyramidBMType::SymmetryType) symmetryArg.getValue());
        matcher->SetAgregator((PyramidBMType::Agregator) agregatorArg.getValue());
        matcher->SetOutputTransformType(PyramidBMType::outRigid);
        matcher->SetAgregThreshold(agregThresholdArg.getValue());
        matcher->SetSeStoppingThreshold(seStoppingThresholdArg.getValue());
        matcher->SetNumberOfPyramidLevels(numPyramidLevelsArg.getValue());
        matcher->SetLastPyramidLevel(lastPyramidLevelArg.getValue());
        matcher->SetVerbose(false);

        if (numThreadsArg.getValue() != 0)
            matcher->SetNumberOfWorkUnits( numThreadsArg.getValue() );

        matcher->SetPercentageKept( percentageKeptArg.getValue() );
        matcher->SetTransformInitializationType(PyramidBMType::GravityCenters);

        matcher->SetFloatingImage(referenceExtractFilter->GetOutput());
        matcher->SetReferenceImage(extractFilter->GetOutput());

        AffineTransformPointer rigidTrsf = AffineTransformType::New();
        rigidTrsf->SetIdentity();
        matcher->SetOutputTransform(rigidTrsf.GetPointer());

        try
        {
            matcher->Update();
        }
        catch (itk::ExceptionObject &e)
        {
            std::cerr << e << std::endl;
            return EXIT_FAILURE;
        }

        rigidTrsf = dynamic_cast <AffineTransformType *> (matcher->GetOutputTransform().GetPointer());

        InputSubImageType::Pointer rigidReference = matcher->GetOutputImage();

        // Then perform directional affine registration
        matcher->SetReferenceImage(rigidReference);
        matcher->SetFloatingImage(extractFilter->GetOutput());
        matcher->SetTransform(PyramidBMType::Directional_Affine);
        matcher->SetOutputTransformType(PyramidBMType::outAffine);
        matcher->SetAffineDirection(directionArg.getValue());
        matcher->SetTransformInitializationType(PyramidBMType::Identity);

        AffineTransformPointer initTrsf = AffineTransformType::New();
        initTrsf->SetIdentity();
        matcher->SetInitialTransform(initTrsf);

        AffineTransformPointer tmpTrsfDirectional = AffineTransformType::New();
        tmpTrsfDirectional->SetIdentity();
        matcher->SetOutputTransform(tmpTrsfDirectional.GetPointer());

        try
        {
            matcher->Update();
        }
        catch (itk::ExceptionObject &e)
        {
            std::cerr << e << std::endl;
            return EXIT_FAILURE;
        }

        // Finally, perform non linear registration to get rid of non linear distortions
        typedef anima::PyramidalDenseSVFMatchingBridge <Dimension> NonLinearPyramidBMType;
        NonLinearPyramidBMType::Pointer nonLinearMatcher = NonLinearPyramidBMType::New();

        nonLinearMatcher->SetReferenceImage(rigidReference);
        nonLinearMatcher->SetFloatingImage(matcher->GetOutputImage());

        // Setting matcher arguments
        nonLinearMatcher->SetBlockSize(blockSizeArg.getValue());
        nonLinearMatcher->SetBlockSpacing(nlBlockSpacingArg.getValue());
        nonLinearMatcher->SetStDevThreshold(stdevThresholdArg.getValue());
        nonLinearMatcher->SetTransform(NonLinearPyramidBMType::Directional_Affine);
        nonLinearMatcher->SetAffineDirection(directionArg.getValue());
        nonLinearMatcher->SetMetric((NonLinearPyramidBMType::Metric) blockMetricArg.getValue());
        nonLinearMatcher->SetOptimizer((NonLinearPyramidBMType::Optimizer) optimizerArg.getValue());
        nonLinearMatcher->SetMaximumIterations(maxIterationsArg.getValue());
        nonLinearMatcher->SetMinimalTransformError(minErrorArg.getValue());
        nonLinearMatcher->SetOptimizerMaximumIterations(optimizerMaxIterationsArg.getValue());
        nonLinearMatcher->SetStepSize(searchStepArg.getValue());
        nonLinearMatcher->SetTranslateUpperBound(translateUpperBoundArg.getValue());
        nonLinearMatcher->SetSymmetryType((NonLinearPyramidBMType::SymmetryType) symmetryArg.getValue());
        nonLinearMatcher->SetAgregator(NonLinearPyramidBMType::Baloo);
        nonLinearMatcher->SetBCHCompositionOrder(1);
        nonLinearMatcher->SetExponentiationOrder(0);
        nonLinearMatcher->SetExtrapolationSigma(extrapolationSigmaArg.getValue());
        nonLinearMatcher->SetElasticSigma(elasticSigmaArg.getValue());
        nonLinearMatcher->SetOutlierSigma(outlierSigmaArg.getValue());
        nonLinearMatcher->SetNumberOfPyramidLevels(numPyramidLevelsArg.getValue());
        nonLinearMatcher->SetLastPyramidLevel(lastPyramidLevelArg.getValue());
        nonLinearMatcher->SetVerbose(false);

        if (numThreadsArg.getValue() != 0)
            nonLinearMatcher->SetNumberOfWorkUnits(numThreadsArg.getValue());

        nonLinearMatcher->SetPercentageKept(percentageKeptArg.getValue());

        try
        {
            nonLinearMatcher->Update();
        }
        catch (itk::ExceptionObject &e)
        {
            std::cerr << e << std::endl;
            return EXIT_FAILURE;
        }

        // Finally, apply transform serie to image
        typedef itk::CompositeTransform <AgregatorType::ScalarType,Dimension> GeneralTransformType;
        GeneralTransformType::Pointer transformSerie = GeneralTransformType::New();
        transformSerie->AddTransform(tmpTrsfDirectional);

        typedef itk::StationaryVelocityFieldTransform <AgregatorType::ScalarType,Dimension> SVFTransformType;
        typedef SVFTransformType::Pointer SVFTransformPointer;

        typedef rpi::DisplacementFieldTransform <AgregatorType::ScalarType,Dimension> DenseTransformType;
        typedef DenseTransformType::Pointer DenseTransformPointer;

        SVFTransformPointer svfPointer = nonLinearMatcher->GetOutputTransform();

        DenseTransformPointer dispTrsf = DenseTransformType::New();
        anima::GetSVFExponential(svfPointer.GetPointer(),dispTrsf.GetPointer(),0,numThreadsArg.getValue(),1.0);

        transformSerie->AddTransform(dispTrsf.GetPointer());

        // Apply rigid matrix to gradient vectors
        AffineTransformType::MatrixType rigidMatrix = rigidTrsf->GetMatrix();
        vnl_vector_fixed <double,3> tmpDir(0.0);
        for (unsigned int j = 0;j < 3;++j)
        {
            for (unsigned int k = 0;k < 3;++k)
                tmpDir[j] += rigidMatrix(j,k) * directions[i][k];
        }

        directions[i] = tmpDir;

        AffineTransformPointer rigidTrsfInverse = AffineTransformType::New();
        rigidTrsf->GetInverse(rigidTrsfInverse);
        transformSerie->AddTransform(rigidTrsfInverse.GetPointer());

        typedef anima::ResampleImageFilter<InputSubImageType, InputSubImageType> ResampleFilterType;
        ResampleFilterType::Pointer scalarResampler = ResampleFilterType::New();

        InputSubImageType::SizeType size = referenceExtractFilter->GetOutput()->GetLargestPossibleRegion().GetSize();
        InputSubImageType::PointType origin = referenceExtractFilter->GetOutput()->GetOrigin();
        InputSubImageType::SpacingType spacing = referenceExtractFilter->GetOutput()->GetSpacing();
        InputSubImageType::DirectionType direction = referenceExtractFilter->GetOutput()->GetDirection();

        scalarResampler->SetTransform(transformSerie);
        scalarResampler->SetSize(size);
        scalarResampler->SetOutputOrigin(origin);
        scalarResampler->SetOutputSpacing(spacing);
        scalarResampler->SetOutputDirection(direction);

        scalarResampler->SetInput(extractFilter->GetOutput());
        if (numThreadsArg.getValue() != 0)
            scalarResampler->SetNumberOfWorkUnits(numThreadsArg.getValue());
        scalarResampler->Update();

        InputSubImageType::RegionType regionSubImage = scalarResampler->GetOutput()->GetLargestPossibleRegion();
        InputImageType::RegionType regionImage = inputImage->GetLargestPossibleRegion();
        regionImage.SetIndex(Dimension,i);
        regionImage.SetSize(Dimension,1);

        InputImageIteratorType outIterator(inputImage,regionImage);
        InputSubImageIteratorType inIterator(scalarResampler->GetOutput(),regionSubImage);

        while (!inIterator.IsAtEnd())
        {
            outIterator.Set(inIterator.Get());

            ++inIterator;
            ++outIterator;
        }
    }

    std::cout << std::endl;

    anima::writeImage <InputImageType> (outArg.getValue(),inputImage);

    // Writing output gradients
    std::ofstream outputFile(outBVecArg.getValue());
    for (unsigned int i = 0;i < 3;++i)
    {
        for (unsigned int j = 0;j < directions.size();++j)
            outputFile << directions[j][i] << " ";

        outputFile << std::endl;
    }

    outputFile.close();

    return EXIT_SUCCESS;
}
