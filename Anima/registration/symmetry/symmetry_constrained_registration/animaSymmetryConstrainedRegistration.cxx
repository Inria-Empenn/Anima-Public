#include <animaPyramidalSymmetryConstrainedRegistrationBridge.h>

#include <itkTimeProbe.h>
#include <animaReadWriteFunctions.h>
#include <itkTransformFileReader.h>

#include <tclap/CmdLine.h>

int main(int ac, const char** av)
{
    // Parsing arguments
    TCLAP::CmdLine  cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    // Setting up parameters
    TCLAP::ValueArg<std::string> fixedArg("r","refimage","Fixed image",true,"","fixed image",cmd);
    TCLAP::ValueArg<std::string> movingArg("m","movingimage","Moving image",true,"","moving image",cmd);

    TCLAP::ValueArg<std::string> fixedSymmetryArg("","ref-sym","Fixed symmetry",true,"","fixed symmetry",cmd);
    TCLAP::ValueArg<std::string> movingSymmetryArg("","moving-sym","Moving symmetry",true,"","moving symmetry",cmd);

    TCLAP::ValueArg<std::string> outArg("o","outputimage","Output (registered) image",true,"","output image",cmd);
    TCLAP::ValueArg<std::string> outputTransformArg("O","outtransform","Output transformation",false,"","output transform",cmd);

    TCLAP::ValueArg<unsigned int> metricArg("","metric","Similarity metric (0: mutual information, 1: normalized mutual information, 2: mean squares, default: 2)",false,2,"similarity metric",cmd);
    TCLAP::SwitchArg fastRegArg("F","fast-reg","registration on the central 2D slice only",cmd,false);

    TCLAP::ValueArg<unsigned int> optimizerMaxIterationsArg("","oi","Maximum iterations for local optimizer (default: 100)",false,100,"maximum local optimizer iterations",cmd);
    TCLAP::ValueArg<unsigned int> histoSizeArg("","hs","Histogram size for mutual information (default: 128)",false,128,"histogram size",cmd);

    TCLAP::ValueArg<double> translateUpperBoundArg("","tub","Upper bound on translation for bobyqa (in voxels, default: 10)",false,10,"Bobyqa translate upper bound",cmd);
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

    typedef anima::PyramidalSymmetryConstrainedRegistrationBridge <double> BridgeType;
    typedef itk::Image <double, 3> InputImageType;

    BridgeType::Pointer matcher = BridgeType::New();

    matcher->SetReferenceImage(anima::readImage <InputImageType> (fixedArg.getValue()));
    matcher->SetdoubleingImage(anima::readImage <InputImageType> (movingArg.getValue()));

    matcher->SetResultFile( outArg.getValue() );
    matcher->SetOutputTransformFile( outputTransformArg.getValue() );

    matcher->SetMetric( (Metric) metricArg.getValue() );
    matcher->SetOptimizerMaximumIterations( optimizerMaxIterationsArg.getValue() );
    matcher->SetUpperBoundAngle(angleUpperBoundArg.getValue() * M_PI / 180.0);
    matcher->SetTranslateUpperBound(translateUpperBoundArg.getValue());
    matcher->SetHistogramSize(histoSizeArg.getValue());
    matcher->SetNumberOfPyramidLevels( numPyramidLevelsArg.getValue() );
    matcher->SetFastRegistration(fastRegArg.isSet());

    if (numThreadsArg.getValue() != 0)
        matcher->SetNumberOfWorkUnits( numThreadsArg.getValue() );

    typedef BridgeType::BaseTransformType BaseTransformType;
    itk::TransformFileReader::Pointer tmpTrRead = itk::TransformFileReader::New();
    tmpTrRead->SetFileName(fixedSymmetryArg.getValue());

    try
    {
        tmpTrRead->Update();

        itk::TransformFileReader::TransformListType trsfList = *(tmpTrRead->GetTransformList());
        itk::TransformFileReader::TransformListType::iterator tr_it = trsfList.begin();

        matcher->SetRefSymmetryTransform (dynamic_cast< BaseTransformType *> ((*tr_it).GetPointer()));
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << "Unable to read reference symmetry transform... " << e << std::endl;
        return EXIT_FAILURE;
    }

    tmpTrRead = itk::TransformFileReader::New();
    tmpTrRead->SetFileName(movingSymmetryArg.getValue());

    try
    {
        tmpTrRead->Update();

        itk::TransformFileReader::TransformListType trsfList = *(tmpTrRead->GetTransformList());
        itk::TransformFileReader::TransformListType::iterator tr_it = trsfList.begin();

        matcher->SetFloSymmetryTransform (dynamic_cast< BaseTransformType *> ((*tr_it).GetPointer()));
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << "Unable to read doubleing symmetry transform... " << e << std::endl;
        return EXIT_FAILURE;
    }

    itk::TimeProbe timer;

    timer.Start();

    try
    {
        matcher->Update();
    }
    catch (itk::ExceptionObject &err)
    {
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }

    matcher->WriteOutputs();

    timer.Stop();

    std::cout << "Elapsed Time: " << timer.GetTotal()  << timer.GetUnit() << std::endl;

    return EXIT_SUCCESS;
}
