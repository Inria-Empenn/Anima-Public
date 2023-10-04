#include <animaMCMTractographyImageFilter.h>
#include <itkTimeProbe.h>

#include <animaMCMFileReader.h>
#include <animaReadWriteFunctions.h>
#include <itkCommand.h>

#include <tclap/CmdLine.h>

#include <animaShapesWriter.h>

//Update progression of the process
void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
    std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
}

int main(int argc,  char*  argv[])
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> mcmArg("i","input-mcm","Input MCM image",true,"","MCM image",cmd);
    TCLAP::ValueArg<std::string> seedMaskArg("s","seed-mask","Seed mask",true,"","seed",cmd);
    TCLAP::ValueArg<std::string> fibersArg("o","fibers","Output fibers",true,"","fibers",cmd);

    TCLAP::ValueArg<std::string> filterMaskArg("","filter-mask","Mask for filtering fibers (default: none)",false,"","filter mask",cmd);
    TCLAP::ValueArg<std::string> cutMaskArg("c","cut-mask","Mask for cutting fibers (default: none)",false,"","cut mask",cmd);
    TCLAP::ValueArg<std::string> forbiddenMaskArg("f","forbidden-mask","Mask for removing fibers (default: none)",false,"","remove mask",cmd);

    TCLAP::ValueArg<double> isoThrArg("I","iso-thr","Isotropic weight threshold (default: 0.8)",false,0.8,"Isotropic weight threshold",cmd);
    TCLAP::ValueArg<double> relWeightArg("R","start-rel-weight","Minimal relative weight of direction to start a fiber (default: 0.2)",false,0.2,"Minimal relative weight of start direction",cmd);

    TCLAP::ValueArg<double> stopAngleArg("a","angle-max","Maximum angle for tracking (default: 60)",false,60.0,"maximum track angle",cmd);
    TCLAP::ValueArg<double> stepLengthArg("","step-length","Length of each step (in mm, default: 0.5)",false,0.5,"step length",cmd);
    TCLAP::ValueArg<int> nbFibersArg("","nb-fibers","Number of starting filters (n*n*n) per voxel (default: 1)",false,1,"number of seeds per voxel",cmd);

    TCLAP::ValueArg<double> minLengthArg("","min-length","Minimum length for a fiber to be considered for computation (default: 10mm)",false,10.0,"minimum length",cmd);
    TCLAP::ValueArg<double> maxLengthArg("","max-length","Maximum length of a tract (default: 200mm)",false,200.0,"maximum length",cmd);

    TCLAP::SwitchArg addLocalDataArg("L","local-data","Add local data information to output tracks",cmd);

    TCLAP::ValueArg<unsigned int> nbThreadsArg("T","nb-threads","Number of threads to run on (default: all available)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef anima::MCMTractographyImageFilter MainFilterType;
    typedef MainFilterType::MaskImageType MaskImageType;
    typedef itk::ImageFileReader <MaskImageType> MaskReaderType;

    MainFilterType::Pointer mcmTracker = MainFilterType::New();

    anima::MCMFileReader <double, 3> mcmReader;
    mcmReader.SetFileName(mcmArg.getValue());
    mcmReader.Update();

    mcmTracker->SetInputImage(mcmReader.GetModelVectorImage());

    mcmTracker->SetSeedingMask(anima::readImage <MaskImageType> (seedMaskArg.getValue()));

    if (filterMaskArg.getValue() != "")
        mcmTracker->SetFilteringMask(anima::readImage <MaskImageType> (filterMaskArg.getValue()));

    if (cutMaskArg.getValue() != "")
        mcmTracker->SetCutMask(anima::readImage <MaskImageType> (cutMaskArg.getValue()));

    if (forbiddenMaskArg.getValue() != "")
        mcmTracker->SetForbiddenMask(anima::readImage <MaskImageType> (forbiddenMaskArg.getValue()));

    mcmTracker->SetNumberOfWorkUnits(nbThreadsArg.getValue());
    mcmTracker->SetNumberOfFibersPerPixel(nbFibersArg.getValue());
    mcmTracker->SetStepProgression(stepLengthArg.getValue());
    mcmTracker->SetStopIsoWeightThreshold(isoThrArg.getValue());
    mcmTracker->SetMinimalDirectionRelativeWeight(relWeightArg.getValue());
    mcmTracker->SetMaxFiberAngle(stopAngleArg.getValue());
    mcmTracker->SetMinLengthFiber(minLengthArg.getValue());
    mcmTracker->SetMaxLengthFiber(maxLengthArg.getValue());

    bool computeLocalColors = (fibersArg.getValue().find(".fds") != std::string::npos) && (addLocalDataArg.isSet());
    mcmTracker->SetComputeLocalColors(computeLocalColors);

    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);
    mcmTracker->AddObserver(itk::ProgressEvent(), callback);

    itk::TimeProbe tmpTime;
    tmpTime.Start();

    mcmTracker->Update();
    mcmTracker->RemoveAllObservers();

    tmpTime.Stop();
    std::cout << "Tracking time: " << tmpTime.GetTotal() << "s" << std::endl;

    anima::ShapesWriter writer;
    writer.SetInputData(mcmTracker->GetOutput());
    writer.SetFileName(fibersArg.getValue());
    writer.Update();

    return EXIT_SUCCESS;
}
