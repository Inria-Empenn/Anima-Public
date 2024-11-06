#include <animaPyramidalSymmetryBridge.h>
#include <animaReadWriteFunctions.h>

#include <itkTimeProbe.h>
#include <tclap/CmdLine.h>

int main(int ac, const char** av)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inputArg("i","ref-image","image",true,"","image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","out-image","Output (registered) image",true,"","output image",cmd);

    TCLAP::ValueArg<std::string> outputRealignTransformArg("","out-realign-trsf","Output transformation to realign on symmetry plane",false,"","output realign transform",cmd);
    TCLAP::ValueArg<std::string> outputTransformArg("O","out-trsf","Output transformation",false,"","output transform",cmd);

    TCLAP::ValueArg<unsigned int> metricArg("","metric","Similarity metric (0: mutual information, 1: mean squares, default: 1)",false,1,"similarity metric",cmd);

    TCLAP::ValueArg<unsigned int> optimizerMaxIterationsArg("","oi","Maximum iterations for local optimizer (default: 100)",false,100,"maximum local optimizer iterations",cmd);
    TCLAP::ValueArg<unsigned int> histoSizeArg("","hs","Histogram size for mutual information (default: 128)",false,128,"histogram size",cmd);

    TCLAP::ValueArg<double> translateUpperBoundArg("","tub","Upper bound on translation for bobyqa (in voxels, default: 6)",false,6,"Bobyqa translate upper bound",cmd);
    TCLAP::ValueArg<double> angleUpperBoundArg("","aub","Upper bound on angles for bobyqa (in degrees, default: 180)",false,180,"Bobyqa angle upper bound",cmd);

    TCLAP::ValueArg<unsigned int> numPyramidLevelsArg("p","pyr","Number of pyramid levels (default: 3)",false,3,"number of pyramid levels",cmd);
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

    typedef anima::PyramidalSymmetryBridge <double, double> PyramidSymType;
    typedef itk::Image <double,3> InputImageType;

    PyramidSymType::Pointer matcher = PyramidSymType::New();

    matcher->SetReferenceImage(anima::readImage <InputImageType> (inputArg.getValue()));
    matcher->SetFloatingImage(anima::readImage <InputImageType> (inputArg.getValue()));

    // set parameters
    matcher->SetMetric((Metric)metricArg.getValue());
    matcher->SetOptimizerMaxIterations(optimizerMaxIterationsArg.getValue());
    matcher->SetHistogramSize(histoSizeArg.getValue());
    matcher->SetUpperBoundDistance(translateUpperBoundArg.getValue());
    matcher->SetUpperBoundAngle(angleUpperBoundArg.getValue() * M_PI / 180.0);
    matcher->SetNumberOfPyramidLevels(numPyramidLevelsArg.getValue());

    if (numThreadsArg.getValue() != 0)
        matcher->SetNumberOfWorkUnits(numThreadsArg.getValue());

    matcher->SetResultFile(outArg.getValue());
    matcher->SetOutputRealignTransformFile(outputRealignTransformArg.getValue());
    matcher->SetOutputTransformFile(outputTransformArg.getValue());

    itk::TimeProbe timer;

    timer.Start();

    matcher->Update();
    matcher->WriteOutputs();

    timer.Stop();

    std::cout << "Elapsed Time: " << timer.GetTotal()  << timer.GetUnit() << std::endl;

    return EXIT_SUCCESS;
}
