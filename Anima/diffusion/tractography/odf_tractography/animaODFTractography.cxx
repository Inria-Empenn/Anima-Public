#include <animaODFTractographyImageFilter.h>
#include <animaGradientFileReader.h>
#include <itkTimeProbe.h>

#include <tclap/CmdLine.h>

#include <animaReadWriteFunctions.h>
#include <animaLogTensorImageFilter.h>

#include <itkCommand.h>

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

    // Mandatory arguments
    TCLAP::ValueArg<std::string> odfArg("i","odf","Input odf image",true,"","odf image",cmd);
    TCLAP::ValueArg<std::string> seedMaskArg("s","seed-mask","Seed mask",true,"","seed",cmd);
    TCLAP::ValueArg<std::string> fibersArg("o","fibers","Output fibers",true,"","fibers",cmd);

    // Optional mask arguments
    TCLAP::ValueArg<std::string> cutMaskArg("c","cut-mask","Mask for cutting fibers (default: none)",false,"","cut mask",cmd);
    TCLAP::ValueArg<std::string> forbiddenMaskArg("f","forbidden-mask","Mask for removing fibers (default: none)",false,"","remove mask",cmd);
    TCLAP::ValueArg<std::string> filterMaskArg("","filter-mask","Mask for filtering fibers (default: none)",false,"","filter mask",cmd);

    TCLAP::ValueArg<double> gfaThrArg("","gfa-thr","GFA threshold (default: 0.2)",false,0.2,"GFA threshold",cmd);
    TCLAP::ValueArg<double> stopAngleArg("a","angle-max","Maximum angle for tracking (default: 60)",false,60.0,"maximum track angle",cmd);
    TCLAP::ValueArg<double> stepLengthArg("","step-length","Length of each step (default: 1)",false,1.0,"step length",cmd);
    //TCLAP::ValueArg<int> nbFibersArg("","nb-fibers","Number of starting particle filters (n*n*n) per voxel (default: 1)",false,1,"number of seeds per voxel",cmd);

    TCLAP::ValueArg<double> minDiffProbaArg("","mdp","Minimal diffusion probability for an ODF direction to be kept (default: 0.02",false,0.01,"Minimal diffusion probability for an ODF direction",cmd);

    TCLAP::ValueArg<double> minLengthArg("","min-length","Minimum length for a fiber to be considered for computation (default: 10mm)",false,10.0,"minimum length",cmd);
    TCLAP::ValueArg<double> maxLengthArg("","max-length","Maximum length of a tract (default: 150mm)",false,150.0,"maximum length",cmd);

    TCLAP::ValueArg<unsigned int> nbThreadsArg("T","nb-threads","Number of threads to run on (default: all available)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    typedef anima::ODFTractographyImageFilter MainFilterType;
    typedef MainFilterType::ModelImageType ModelImageType;
    typedef MainFilterType::MaskImageType MaskImageType;
    typedef MainFilterType::Vector3DType Vector3DType;

    MainFilterType::Pointer odfTracker = MainFilterType::New();

    odfTracker->SetNumberOfWorkUnits(nbThreadsArg.getValue());
    odfTracker->SetInputImage(anima::readImage <ModelImageType> (odfArg.getValue()));

    // Load seed mask
    odfTracker->SetSeedingMask(anima::readImage <MaskImageType> (seedMaskArg.getValue()));

    // Load cut mask
    if (cutMaskArg.getValue() != "")
        odfTracker->SetCutMask(anima::readImage <MaskImageType> (cutMaskArg.getValue()));

    // Load forbidden mask
    if (forbiddenMaskArg.getValue() != "")
        odfTracker->SetForbiddenMask(anima::readImage <MaskImageType> (forbiddenMaskArg.getValue()));

    // Load filter mask
    if (filterMaskArg.getValue() != "")
        odfTracker->SetFilteringMask(anima::readImage <MaskImageType> (filterMaskArg.getValue()));

    odfTracker->SetStepProgression(stepLengthArg.getValue());
    odfTracker->SetGFAThreshold(gfaThrArg.getValue());
    odfTracker->SetMaxFiberAngle(stopAngleArg.getValue());
    odfTracker->SetMinLengthFiber(minLengthArg.getValue());
    odfTracker->SetMaxLengthFiber(maxLengthArg.getValue());
    odfTracker->SetMinimalDiffusionProbability(minDiffProbaArg.getValue());

    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);
    odfTracker->AddObserver(itk::ProgressEvent(), callback);

    itk::TimeProbe tmpTime;
    tmpTime.Start();

    odfTracker->Update();
    odfTracker->RemoveAllObservers();

    tmpTime.Stop();
    std::cout << "Tracking time: " << tmpTime.GetTotal() << "s" << std::endl;

    anima::ShapesWriter writer;
    writer.SetInputData(odfTracker->GetOutput());
    writer.SetFileName(fibersArg.getValue());
    writer.Update();

    return EXIT_SUCCESS;
}
