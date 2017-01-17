#include <animaODFProbabilisticTractographyImageFilter.h>
#include <animaGradientFileReader.h>
#include <itkTimeProbe.h>

#include <tclap/CmdLine.h>

#include <itkImageFileReader.h>
#include <itkCommand.h>

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

    // Mandatory arguments
    TCLAP::ValueArg<std::string> dwiArg("i","dwi","Input 4D diffusion image",true,"","dwi image",cmd);
    TCLAP::ValueArg<std::string> seedMaskArg("s","seed-mask","Seed mask",true,"","seed",cmd);
    TCLAP::ValueArg<std::string> fibersArg("o","fibers","Output fibers",true,"","fibers",cmd);
    TCLAP::ValueArg<std::string> gradientsArg("g","gradients","Gradient table",true,"","gradients",cmd);
    TCLAP::ValueArg<std::string> bvaluesArg("b","bvalues","B-value list",true,"","b-values",cmd);

    TCLAP::ValueArg<int> colinearityModeArg("","col-init-mode",
                                            "Colinearity mode for initialization - 0: center, 1: outward, 2: top, 3: bottom, 4: left, 5: right, 6: front, 7: back (default: 0)",
                                            false,0,"colinearity mode for initialization",cmd);
    TCLAP::ValueArg<int> initialDirectionModeArg("","init-mode",
                                                 "Mode for initialization - 0: take most colinear direction, 1: take highest weighted direction (default: 1)",
                                                 false,1,"mode for initialization",cmd);

    // Optional mask arguments
    TCLAP::ValueArg<std::string> cutMaskArg("c","cut-mask","Mask for cutting fibers (default: none)",false,"","cut mask",cmd);
    TCLAP::ValueArg<std::string> forbiddenMaskArg("f","forbidden-mask","Mask for removing fibers (default: none)",false,"","remove mask",cmd);
    TCLAP::ValueArg<std::string> filterMaskArg("","filter-mask","Mask for filtering fibers (default: none)",false,"","filter mask",cmd);

    TCLAP::ValueArg<double> gfaThrArg("","gfa-thr","GFA threshold (default: 0.2)",false,0.2,"GFA threshold",cmd);
    TCLAP::ValueArg<double> stepLengthArg("","step-length","Length of each step (default: 1)",false,1.0,"step length",cmd);
    TCLAP::ValueArg<int> nbFibersArg("","nb-fibers","Number of starting particle filters (n*n*n) per voxel (default: 1)",false,1,"number of seeds per voxel",cmd);

    TCLAP::ValueArg<double> minLengthArg("","min-length","Minimum length for a fiber to be considered for computation (default: 10mm)",false,10.0,"minimum length",cmd);
    TCLAP::ValueArg<double> maxLengthArg("","max-length","Maximum length of a tract (default: 150mm)",false,150.0,"maximum length",cmd);

    TCLAP::ValueArg<unsigned int> nbParticlesArg("n","nb-particles","Number of particles per filter (default: 1000)",false,1000,"number of particles",cmd);
    TCLAP::ValueArg<unsigned int> clusterMinSizeArg("","cluster-min-size","Minimal number of particles per cluster before split (default: 10)",false,10,"minimal cluster size",cmd);

    TCLAP::ValueArg<double> resampThrArg("r","resamp-thr","Resampling relative threshold (default: 0.8)",false,0.8,"resampling threshold",cmd);

    TCLAP::ValueArg<double> minDiffProbaArg("","mdp","Minimal diffusion probability for an ODF direction to be kept (default: 0.02",false,0.01,"Minimal diffusion probability for an ODF direction",cmd);
    TCLAP::ValueArg<double> trashThrArg("","trash-thr","Relative threshold to keep fibers in trash (default: 0.1)",false,0.1,"trash threshold",cmd);
    TCLAP::ValueArg<double> kappaPriorArg("k","kappa-prior","Kappa of prior distribution (default: 15)",false,15.0,"prior kappa",cmd);
    TCLAP::ValueArg<unsigned int> odfOrderArg("K","odf-order","ODF order (default: 4)",false,4,"ODF order",cmd);
    TCLAP::ValueArg<double> curvScaleArg("","cs","Scale for ODF curvature to get vMF kappa (default: 6)",false,6.0,"scale for ODF curvature",cmd);

    TCLAP::ValueArg<double> distThrArg("","dist-thr","Hausdorff distance threshold for mergine clusters (default: 0.5)",false,0.5,"merging threshold",cmd);
    TCLAP::ValueArg<double> kappaThrArg("","kappa-thr","Kappa threshold for splitting clusters (default: 30)",false,30.0,"splitting threshold",cmd);

    TCLAP::ValueArg<unsigned int> clusterDistArg("","cluster-dist","Distance between clusters: choices are 0 (AHD), 1 (HD, default) or 2 (MHD)",false,1,"cluster distance",cmd);

    TCLAP::SwitchArg averageClustersArg("M","average-clusters","Output only cluster mean",cmd,false);

    TCLAP::ValueArg<unsigned int> nbThreadsArg("t","nb-threads","Number of threads to run on (default: all available)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }
    
    typedef anima::ODFProbabilisticTractographyImageFilter MainFilterType;
    typedef MainFilterType::Input4DImageType Input4DImageType;
    typedef itk::ImageFileReader <Input4DImageType> Input4DReaderType;
    typedef MainFilterType::MaskImageType MaskImageType;
    typedef itk::ImageFileReader <MaskImageType> MaskReaderType;
    typedef MainFilterType::Vector3DType Vector3DType;
    
    MainFilterType::Pointer odfTracker = MainFilterType::New();

    Input4DReaderType::Pointer rawReader = Input4DReaderType::New();
    rawReader->SetFileName(dwiArg.getValue());
    std::cout << "Loading the diffusion images... " << std::endl;
    rawReader->Update();
    odfTracker->SetNumberOfThreads(nbThreadsArg.getValue());
    odfTracker->SetInputImagesFrom4DImage(rawReader->GetOutput());
    rawReader = 0;

    // Read gradient table and b-value list
    typedef anima::GradientFileReader < Vector3DType, double > GFReaderType;
    GFReaderType gfReader;
    gfReader.SetGradientFileName(gradientsArg.getValue());
    gfReader.SetBValueBaseString(bvaluesArg.getValue());
    gfReader.SetGradientIndependentNormalization(false);
    
    gfReader.Update();
    
    GFReaderType::GradientVectorType directions = gfReader.GetGradients();
    
    for(unsigned int i = 0;i < directions.size();++i)
        odfTracker->AddGradientDirection(i,directions[i]);
    
    GFReaderType::BValueVectorType mb = gfReader.GetBValues();
    
    odfTracker->SetBValuesList(mb);

    odfTracker->SetInitialColinearityDirection((MainFilterType::ColinearityDirectionType)colinearityModeArg.getValue());
    odfTracker->SetInitialDirectionMode((MainFilterType::InitialDirectionModeType)initialDirectionModeArg.getValue());
    
    MaskReaderType::Pointer roiReader = MaskReaderType::New();
    roiReader->SetFileName(seedMaskArg.getValue());
    roiReader->Update();
    
    odfTracker->SetSeedMask(roiReader->GetOutput());
    
    if (cutMaskArg.getValue() != "")
    {
        MaskReaderType::Pointer cutReader = MaskReaderType::New();
        cutReader->SetFileName(cutMaskArg.getValue());
        cutReader->Update();
        
        odfTracker->SetCutMask(cutReader->GetOutput());
    }
    
    if (forbiddenMaskArg.getValue() != "")
    {
        MaskReaderType::Pointer forbiddenReader = MaskReaderType::New();
        forbiddenReader->SetFileName(forbiddenMaskArg.getValue());
        forbiddenReader->Update();
        
        odfTracker->SetForbiddenMask(forbiddenReader->GetOutput());
    }
    
    if (filterMaskArg.getValue() != "")
    {
        MaskReaderType::Pointer filterMaskReader = MaskReaderType::New();
        filterMaskReader->SetFileName(filterMaskArg.getValue());
        filterMaskReader->Update();
        
        odfTracker->SetFilterMask(filterMaskReader->GetOutput());
    }

    odfTracker->SetNumberOfFibersPerPixel(nbFibersArg.getValue());
    odfTracker->SetStepProgression(stepLengthArg.getValue());
    odfTracker->SetGFAThreshold(gfaThrArg.getValue());
    odfTracker->SetMinLengthFiber(minLengthArg.getValue());
    odfTracker->SetMaxLengthFiber(maxLengthArg.getValue());
    
    odfTracker->SetNumberOfParticles(nbParticlesArg.getValue());
    odfTracker->SetMinimalNumberOfParticlesPerClass(clusterMinSizeArg.getValue());
    odfTracker->SetResamplingThreshold(resampThrArg.getValue());
    odfTracker->SetFiberTrashThreshold(trashThrArg.getValue());
    
    odfTracker->SetMinimalDiffusionProbability(minDiffProbaArg.getValue());
    odfTracker->SetKappaOfPriorDistribution(kappaPriorArg.getValue());
    odfTracker->SetODFSHOrder(odfOrderArg.getValue());
    
    odfTracker->SetPositionDistanceFuseThreshold(distThrArg.getValue());
    odfTracker->SetKappaSplitThreshold(kappaThrArg.getValue());
    odfTracker->SetClusterDistance(clusterDistArg.getValue());
    odfTracker->SetCurvatureScale(curvScaleArg.getValue());
    
    odfTracker->SetComputeLocalColors(fibersArg.getValue().find(".fds") == std::string::npos);
    odfTracker->SetMAPMergeFibers(averageClustersArg.getValue());
    
    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);
    odfTracker->AddObserver(itk::ProgressEvent(), callback);

    itk::TimeProbe tmpTime;
    tmpTime.Start();
    
    odfTracker->Update();
    odfTracker->RemoveAllObservers();

    tmpTime.Stop();
    std::cout << "Tracking time: " << tmpTime.GetTotal() << "s" << std::endl;
    
    anima::FibersWriter writer;
    writer.SetInputData(odfTracker->GetOutput());
    writer.SetFileName(fibersArg.getValue());
    writer.Update();
    
    return 0;
}
