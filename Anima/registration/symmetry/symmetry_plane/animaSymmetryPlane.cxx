#include <animaPyramidalSymmetryBridge.h>

#include <itkTimeProbe.h>
#include <tclap/CmdLine.h>

int main(int ac, const char** av)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inputArg("i","ref-image","image",true,"","image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","out-image","Output (registered) image",true,"","output image",cmd);

    TCLAP::ValueArg<std::string> outputRealignTransformArg("","out-realign-trsf","Output transformation to realign on symmetry plane",false,"","output realign transform",cmd);
    TCLAP::ValueArg<std::string> outputTransformArg("O","out-trsf","Output transformation",false,"","output transform",cmd);

    TCLAP::ValueArg<unsigned int> metricArg("","metric","Similarity metric (0: mutual information, 1: mean squares, default: 1)",false,1,"similarity metric",cmd);
    TCLAP::ValueArg<unsigned int> optimizerArg("","opt","Optimizer for trsf search (0: Newuoa, 1: Powell, default: 0)",false,0,"optimizer",cmd);

    TCLAP::ValueArg<unsigned int> optimizerMaxIterationsArg("","oi","Maximum iterations for local optimizer (default: 100)",false,100,"maximum local optimizer iterations",cmd);
    TCLAP::ValueArg<unsigned int> histoSizeArg("","hs","Histogram size for mutual information (default: 128)",false,128,"histogram size",cmd);

    TCLAP::ValueArg<double> searchRadiusArg("","sr","Search radius in pixels (exhaustive search window, rho start for bobyqa, default: 2)",false,2,"optimizer search radius",cmd);
    TCLAP::ValueArg<double> searchAngleRadiusArg("","sar","Search angle radius in degrees (rho start for bobyqa, default: 5)",false,5,"optimizer search angle radius",cmd);
    TCLAP::ValueArg<double> finalRadiusArg("","fr","Final radius (rho end for bobyqa, default: 0.001)",false,0.001,"optimizer final radius",cmd);

    TCLAP::ValueArg<unsigned int> numPyramidLevelsArg("p","pyr","Number of pyramid levels (default: 3)",false,3,"number of pyramid levels",cmd);
    TCLAP::ValueArg<unsigned int> numThreadsArg("T","threads","Number of execution threads (default: 0 = all cores)",false,0,"number of threads",cmd);

    try
    {
        cmd.parse(ac,av);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    typedef anima::PyramidalSymmetryBridge <float, double> PyramidSymType;
    typedef itk::Image <float,3> InputImageType;

    PyramidSymType::Pointer matcher = PyramidSymType::New();
    typedef itk::ImageFileReader<InputImageType> ReaderType;

    ReaderType::Pointer tmpReadRef = ReaderType::New();
    tmpReadRef->SetFileName(inputArg.getValue());
    tmpReadRef->Update();

    matcher->SetReferenceImage(tmpReadRef->GetOutput());

    // Cheating here, the two images ar the same
    ReaderType::Pointer tmpReadFlo = ReaderType::New();
    tmpReadFlo->SetFileName(inputArg.getValue());
    tmpReadFlo->Update();

    matcher->SetFloatingImage(tmpReadFlo->GetOutput());

    // set parameters
    matcher->SetMetric((Metric)metricArg.getValue());
    matcher->SetOptimizerType((OptimizerType)optimizerArg.getValue());
    matcher->SetOptimizerMaxIterations(optimizerMaxIterationsArg.getValue());
    matcher->SetHistogramSize(histoSizeArg.getValue());
    matcher->SetSearchRadius(searchRadiusArg.getValue());
    matcher->SetSearchAngleRadius(searchAngleRadiusArg.getValue());
    matcher->SetFinalRadius(finalRadiusArg.getValue());
    matcher->SetNumberOfPyramidLevels(numPyramidLevelsArg.getValue());

    if (numThreadsArg.getValue() != 0)
        matcher->SetNumberOfThreads(numThreadsArg.getValue());

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
