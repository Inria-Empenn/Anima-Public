#include <animaMCMProbabilisticTractographyImageFilter.h>
#include <animaGradientFileReader.h>
#include <itkTimeProbe.h>

#include <tclap/CmdLine.h>
#include <itkCommand.h>

#include <animaReadWriteFunctions.h>
#include <animaMCMFileReader.h>

#include <animaShapesWriter.h>

void ComputeKappaPolynomialCoefficients(std::vector <double> &resVal)
{
    // Values fitted from [Zhang, 2009, Fig.3] : NEEDS TO BE RE-CHECKED
    resVal[2] = 198.81;
    resVal[1] = -73.214;
    resVal[0] = 10.976;
}

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
    TCLAP::ValueArg<std::string> mcmArg("i","mcm","MCM image",true,"","mcm image",cmd);
    TCLAP::ValueArg<std::string> b0Arg("b","b0","B0 image",true,"","B0 image",cmd);
    TCLAP::ValueArg<std::string> noiseArg("N","noise","Noise image",true,"","noise image",cmd);

    TCLAP::ValueArg<std::string> seedMaskArg("s","seed-mask","Seed mask",true,"","seed",cmd);
    TCLAP::ValueArg<std::string> fibersArg("o","fibers","Output fibers",true,"","fibers",cmd);
    
    TCLAP::ValueArg<unsigned int> colinearityModeArg("","col-init-mode",
                                                     "Colinearity mode for initialization - 0: center, 1: outward, 2: top, 3: bottom, 4: left, 5: right, 6: front, 7: back (default: 0)",
                                                     false,0,"colinearity mode for initialization",cmd);
    TCLAP::ValueArg<unsigned int> initialDirectionModeArg("","init-mode",
                                                          "Mode for initialization - 0: take most colinear direction, 1: take highest weighted direction (default: 1)",
                                                          false,1,"mode for initialization",cmd);

    // Optional mask arguments
    TCLAP::ValueArg<std::string> cutMaskArg("c","cut-mask","Mask for cutting fibers (default: none)",false,"","cut mask",cmd);
    TCLAP::ValueArg<std::string> forbiddenMaskArg("f","forbidden-mask","Mask for removing fibers (default: none)",false,"","remove mask",cmd);
    TCLAP::ValueArg<std::string> filterMaskArg("","filter-mask","Mask for filtering fibers (default: none)",false,"","filter mask",cmd);
    
    TCLAP::ValueArg<double> faThrArg("","fa-thr","FA threshold (default: 0.5)",false,0.5,"fa threshold",cmd);
    TCLAP::ValueArg<double> stepLengthArg("","step-length","Length of each step (default: 1)",false,1.0,"step length",cmd);
    TCLAP::ValueArg<int> nbFibersArg("","nb-fibers","Number of starting particle filters (n*n*n) per voxel (default: 1)",false,1,"number of seeds per voxel",cmd);
    
    TCLAP::ValueArg<double> minLengthArg("","min-length","Minimum length for a fiber to be considered for computation (default: 10mm)",false,10.0,"minimum length",cmd);
    TCLAP::ValueArg<double> maxLengthArg("","max-length","Maximum length of a tract (default: 150mm)",false,150.0,"maximum length",cmd);
    
    TCLAP::ValueArg<unsigned int> nbParticlesArg("n","nb-particles","Number of particles per filter (default: 1000)",false,1000,"number of particles",cmd);
    TCLAP::ValueArg<unsigned int> clusterMinSizeArg("","cluster-min-size","Minimal number of particles per cluster before split (default: 10)",false,10,"minimal cluster size",cmd);
    
    TCLAP::ValueArg<double> resampThrArg("r","resamp-thr","Resampling relative threshold (default: 0.8)",false,0.8,"resampling threshold",cmd);
    
    TCLAP::ValueArg<double> trashThrArg("","trash-thr","Relative threshold to keep fibers in trash (default: 0.1)",false,0.1,"trash threshold",cmd);
    TCLAP::ValueArg<double> kappaPriorArg("k","kappa-prior","Kappa of prior distribution (default: 15)",false,15.0,"prior kappa",cmd);

    TCLAP::ValueArg<double> distThrArg("","dist-thr","Hausdorff distance threshold for mergine clusters (default: 0.5)",false,0.5,"merging threshold",cmd);
    TCLAP::ValueArg<double> kappaThrArg("","kappa-thr","Kappa threshold for splitting clusters (default: 30)",false,30.0,"splitting threshold",cmd);
    
    TCLAP::ValueArg<unsigned int> clusterDistArg("","cluster-dist","Distance between clusters: choices are 0 (AHD), 1 (HD, default) or 2 (MHD)",false,1,"cluster distance",cmd);
    
    TCLAP::ValueArg<double> freeWaterFracArg("","free-water-frac","Free water fraction threshold to stop fibers (default: 0.8)",false,0.8,"free water frac thr",cmd);
    
    TCLAP::SwitchArg averageClustersArg("M","average-clusters","Output only cluster mean",cmd,false);
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
    
    typedef anima::MCMProbabilisticTractographyImageFilter MainFilterType;
    typedef MainFilterType::InputModelImageType InputModelImageType;
    typedef MainFilterType::MaskImageType MaskImageType;
    typedef MainFilterType::Vector3DType Vector3DType;
    
    MainFilterType::Pointer mcmTracker = MainFilterType::New();

    anima::MCMFileReader <double,3> mcmReader;
    mcmReader.SetFileName(mcmArg.getValue());
    mcmReader.Update();

    mcmTracker->SetInputModelImage(mcmReader.GetModelVectorImage());

    std::vector <double> kappaCoefficients(3,0);
    ComputeKappaPolynomialCoefficients(kappaCoefficients);
    mcmTracker->SetKappaPolynomialCoefficients(kappaCoefficients);
        
    mcmTracker->SetInitialColinearityDirection((MainFilterType::ColinearityDirectionType)colinearityModeArg.getValue());
    mcmTracker->SetInitialDirectionMode((MainFilterType::InitialDirectionModeType)initialDirectionModeArg.getValue());

    // Load seed mask
    mcmTracker->SetSeedMask(anima::readImage <MaskImageType> (seedMaskArg.getValue()));

    // Load cut mask
    if (cutMaskArg.getValue() != "")
        mcmTracker->SetCutMask(anima::readImage <MaskImageType> (cutMaskArg.getValue()));

    // Load forbidden mask
    if (forbiddenMaskArg.getValue() != "")
        mcmTracker->SetForbiddenMask(anima::readImage <MaskImageType> (forbiddenMaskArg.getValue()));

    // Load filter mask
    if (filterMaskArg.getValue() != "")
        mcmTracker->SetFilterMask(anima::readImage <MaskImageType> (filterMaskArg.getValue()));

    typedef MainFilterType::ScalarImageType ScalarImageType;
    mcmTracker->SetB0Image(anima::readImage <ScalarImageType> (b0Arg.getValue()));
    mcmTracker->SetNoiseImage(anima::readImage <ScalarImageType> (noiseArg.getValue()));

    mcmTracker->SetNumberOfWorkUnits(nbThreadsArg.getValue());
    mcmTracker->SetNumberOfFibersPerPixel(nbFibersArg.getValue());
    mcmTracker->SetStepProgression(stepLengthArg.getValue());
    mcmTracker->SetFAThreshold(faThrArg.getValue());
    mcmTracker->SetMinLengthFiber(minLengthArg.getValue());
    mcmTracker->SetMaxLengthFiber(maxLengthArg.getValue());
    
    mcmTracker->SetNumberOfParticles(nbParticlesArg.getValue());
    mcmTracker->SetMinimalNumberOfParticlesPerClass(clusterMinSizeArg.getValue());
    mcmTracker->SetResamplingThreshold(resampThrArg.getValue());
    mcmTracker->SetFiberTrashThreshold(trashThrArg.getValue());
    
    mcmTracker->SetIsotropicThreshold(freeWaterFracArg.getValue());
    mcmTracker->SetKappaOfPriorDistribution(kappaPriorArg.getValue());

    mcmTracker->SetPositionDistanceFuseThreshold(distThrArg.getValue());
    mcmTracker->SetKappaSplitThreshold(kappaThrArg.getValue());
    mcmTracker->SetClusterDistance(clusterDistArg.getValue());
    
    bool computeLocalColors = (fibersArg.getValue().find(".fds") != std::string::npos) && (addLocalDataArg.isSet());
    mcmTracker->SetComputeLocalColors(computeLocalColors);
    mcmTracker->SetMAPMergeFibers(averageClustersArg.isSet());

    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);
    mcmTracker->AddObserver(itk::ProgressEvent(), callback);

    itk::TimeProbe tmpTime;
    tmpTime.Start();
    
    try
    {
        mcmTracker->Update();
        mcmTracker->RemoveAllObservers();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }
    
    tmpTime.Stop();
    std::cout << "Tracking time: " << tmpTime.GetTotal() << "s" << std::endl;
    
    anima::ShapesWriter writer;
    writer.SetInputData(mcmTracker->GetOutput());
    writer.SetFileName(fibersArg.getValue());
    writer.Update();
    
    return EXIT_SUCCESS;
}
