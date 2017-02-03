#include <animaDTITractographyImageFilter.h>
#include <itkTimeProbe.h>

#include <animaReadWriteFunctions.h>
#include <animaLogTensorImageFilter.h>
#include <itkMultiThreader.h>
#include <itkCommand.h>

#include <tclap/CmdLine.h>

#include <animaFibersWriter.h>

//Update progression of the process
void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
    std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
}

int main(int argc,  char*  argv[])
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> dtiArg("i","dti","Input DT image",true,"","DT image",cmd);
    TCLAP::ValueArg<std::string> seedMaskArg("s","seed-mask","Seed mask",true,"","seed",cmd);
    TCLAP::ValueArg<std::string> fibersArg("o","fibers","Output fibers",true,"","fibers",cmd);

    TCLAP::ValueArg<std::string> cutMaskArg("c","cut-mask","Mask for cutting fibers (default: none)",false,"","cut mask",cmd);
    TCLAP::ValueArg<std::string> forbiddenMaskArg("f","forbidden-mask","Mask for removing fibers (default: none)",false,"","remove mask",cmd);

    TCLAP::ValueArg<double> faThrArg("","fa-thr","FA threshold (default: 0.1)",false,0.1,"fa threshold",cmd);
    TCLAP::ValueArg<double> minNewModelWeightArg("w","min-weight","Minimal model direction weight wrt previous (default: 0.25)",false,0.25,"minimal model direction weight",cmd);

    TCLAP::ValueArg<double> stopAngleArg("a","angle-max","Maximum angle for tracking (default: 60)",false,60.0,"maximum track angle",cmd);
    TCLAP::ValueArg<double> stepLengthArg("","step-length","Length of each step (default: 1)",false,1.0,"step length",cmd);
    TCLAP::ValueArg<int> nbFibersArg("","nb-fibers","Number of starting filters (n*n*n) per voxel (default: 1)",false,1,"number of seeds per voxel",cmd);

    TCLAP::ValueArg<double> minLengthArg("","min-length","Minimum length for a fiber to be considered for computation (default: 10mm)",false,10.0,"minimum length",cmd);
    TCLAP::ValueArg<double> maxLengthArg("","max-length","Maximum length of a tract (default: 150mm)",false,150.0,"maximum length",cmd);

    TCLAP::ValueArg<unsigned int> nbThreadsArg("T","nb-threads","Number of threads to run on (default: all available)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    typedef anima::dtiTractographyImageFilter MainFilterType;
    typedef MainFilterType::ModelImageType DTIImageType;
    typedef MainFilterType::MaskImageType MaskImageType;

    MainFilterType::Pointer dtiTracker = MainFilterType::New();

    typedef anima::LogTensorImageFilter <DTIImageType::IOPixelType,DTIImageType::ImageDimension> LogTensorFilterType;
    LogTensorFilterType::Pointer tensorLogger = LogTensorFilterType::New();

    tensorLogger->SetInput(anima::readImage <DTIImageType> (dtiArg.getValue()));
    tensorLogger->SetScaleNonDiagonal(false);
    tensorLogger->SetNumberOfThreads(nbThreadsArg.getValue());

    tensorLogger->Update();

    dtiTracker->SetInputImage(tensorLogger->GetOutput());

    dtiTracker->SetTrackingMask(anima::readImage <MaskImageType> (seedMaskArg.getValue()));

    if (cutMaskArg.getValue() != "")
        dtiTracker->SetCutMask(anima::readImage <MaskImageType> (cutMaskArg.getValue()));

    if (forbiddenMaskArg.getValue() != "")
        dtiTracker->SetForbiddenMask(anima::readImage <MaskImageType> (forbiddenMaskArg.getValue()));

    dtiTracker->SetNumberOfThreads(nbThreadsArg.getValue());
    dtiTracker->SetNumberOfFibersPerPixel(nbFibersArg.getValue());
    dtiTracker->SetStepProgression(stepLengthArg.getValue());
    dtiTracker->SetStopFAThreshold(faThrArg.getValue());
    dtiTracker->SetMinimalModelWeight(minNewModelWeightArg.getValue());
    dtiTracker->SetMaxFiberAngle(stopAngleArg.getValue());
    dtiTracker->SetMinLengthFiber(minLengthArg.getValue());
    dtiTracker->SetMaxLengthFiber(maxLengthArg.getValue());
    dtiTracker->SetComputeLocalColors(!fibersArg.getValue().find(".fds"));

    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);
    dtiTracker->AddObserver(itk::ProgressEvent(), callback);

    itk::TimeProbe tmpTime;
    tmpTime.Start();

    dtiTracker->Update();
    dtiTracker->RemoveAllObservers();

    tmpTime.Stop();
    std::cout << "Tracking time: " << tmpTime.GetTotal() << "s" << std::endl;

    anima::FibersWriter writer;
    writer.SetInputData(dtiTracker->GetOutput());
    writer.SetFileName(fibersArg.getValue());
    writer.Update();

    return 0;
}
