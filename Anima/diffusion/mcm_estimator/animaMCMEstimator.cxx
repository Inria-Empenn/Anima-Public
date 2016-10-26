#include <fstream>
#include <itkTimeProbe.h>

#include <tclap/CmdLine.h>

#include <animaMCMEstimatorImageFilter.h>
#include <animaGradientFileReader.h>
#include <animaVectorOperations.h>
#include <animaReadWriteFunctions.h>

//Update progression of the process
void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
    std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
}

int main(int argc,  char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',ANIMA_VERSION);

    // Required arguments
    TCLAP::ValueArg<std::string> dwiArg("i", "dwi", "DWI 4D Volume", true, "", "DWI images", cmd);
    TCLAP::ValueArg<std::string> gradsArg("g", "grads", "Gradient table", true, "", "gradients", cmd);
    TCLAP::ValueArg<std::string> bvalsArg("b", "bvals", "B-value list", true, "", "b-values", cmd);
    TCLAP::SwitchArg bvalueScaleArg("B","b-no-scale","Do not scale b-values according to gradient norm",cmd);

    // Outputs
    TCLAP::ValueArg<std::string> outArg("o", "out-mcm", "MCM output volume", true, "", "MCM output", cmd);
    TCLAP::ValueArg<std::string> aicArg("a", "out-aic", "Output estimated AICu image", false, "", "AICu output image", cmd);
    TCLAP::ValueArg<std::string> outB0Arg("", "out-b0", "Output estimated B0 image", false, "", "B0 output image", cmd);
    TCLAP::ValueArg<std::string> outSigmaArg("", "out-sig", "Output estimated noise sigma square image", false, "", "noise sigma square output image", cmd);
    TCLAP::ValueArg<std::string> outMoseArg("", "out-mose", "Output model selection map", false, "", "model selection output image", cmd);
    TCLAP::ValueArg<std::string> inMoseArg("I", "in-mose", "Input model selection map (overrides model selection step)", false, "", "model selection input image", cmd);

    // Optional arguments
    TCLAP::ValueArg<std::string> computationMaskArg("m", "mask", "Computation mask", false, "", "computation mask", cmd);
    TCLAP::ValueArg<double> b0thrArg("", "b0-thr", "Background threshold on B0 value (default: 10)", false, 10.0, "B0 theshold", cmd);

    TCLAP::ValueArg<unsigned int> nbFasciclesArg("n", "nb-fascicles", "Number of computed fascicles (default: 2)", false, 2, "number of fascicles", cmd);
    TCLAP::ValueArg<unsigned int> compartmentTypeArg("c", "comp-type", "Compartment type for fascicles: 1: stick, 2: zeppelin, 3: tensor, 4: noddi (default: 3)", false, 3, "fascicles type", cmd);
    TCLAP::SwitchArg aicSelectNbCompartmentsArg("M", "opt-nb-comp", "Activate AICC-based number of compartments selection", cmd, false);

    TCLAP::SwitchArg vnlDerivArg("D", "vnl-deriv", "Compute derivative by forward difference (used only with levenberg-marquardt algorithm)", cmd, false);
    
    TCLAP::SwitchArg freeWaterCompartmentArg("F", "free-water", "Model with free water", cmd, false);
    TCLAP::SwitchArg stationaryWaterCompartmentArg("S", "stationary-water", "Model with stationary water", cmd, false);
    TCLAP::SwitchArg restrictedWaterCompartmentArg("R", "restricted-water", "Model with restricted water", cmd, false);

    TCLAP::SwitchArg fixWeightsArg("", "fix-weights", "Fix compartment weights", cmd, false);
    TCLAP::ValueArg<double> freeWaterWeightArg("f", "fw-weight", "Free water weight, used if free water compartment (default: 0.1)", false, 0.1, "free water weight", cmd);
    TCLAP::ValueArg<double> stationaryWaterWeightArg("s", "sw-weight", "Stationary water weight, used if stationary water compartment (default: 0.05)", false, 0.05, "stationary water weight", cmd);
    TCLAP::ValueArg<double> restrictedWaterWeightArg("r", "rw-weight", "Restricted water weight, used if restricted water compartment (default: 0.1)", false, 0.1, "restricted water weight", cmd);
    TCLAP::SwitchArg fixDiffArg("", "fix-diff", "Fix diffusivity value", cmd, false);
    TCLAP::SwitchArg optFWDiffArg("", "opt-free-water-diff", "Optimize free water diffusivity value", cmd, false);
    TCLAP::SwitchArg optIRWDiffArg("", "opt-ir-water-diff", "Optimize isotropic restricted water diffusivity value", cmd, false);

    TCLAP::SwitchArg commonDiffusivitiesArg("", "common-diffusivities", "Share diffusivity values among compartments", cmd, false);
    TCLAP::SwitchArg diffInitAbsoluteArg("", "diff-init-abs", "Do not compute initial stick diffusivities from DTI (take absolute values)", cmd, false);

    TCLAP::ValueArg<std::string> optimizerArg("", "optimizer", "Optimizer for estimation: bobyqa (default), ccsaq, bfgs or levenberg", false, "bobyqa", "optimizer", cmd);
    TCLAP::ValueArg<unsigned int> nbRandomRestartsArg("", "random-restarts", "Number of random restarts when searching for more than 2 fascicles (default: 6)", false, 6, "number of random restarts", cmd);
    TCLAP::ValueArg<double> absCostChangeArg("", "abs-cost-change", "Cost function change to stop estimation (default: 0.01)", false, 0.01, "cost change threshold", cmd);
    TCLAP::ValueArg <unsigned int> mlModeArg("", "ml-mode", "ML estimation strategy: marginal likelihood (0), profile likelihood (1, default), Variable projection (2)", false, 1, "ML mode", cmd);
    TCLAP::ValueArg<double> xTolArg("x", "x-tol", "Tolerance for position in optimization (default: 0 -> 1.0e-4 or 1.0e-7 for bobyqa)", false, 0, "position tolerance", cmd);
    TCLAP::ValueArg<double> gTolArg("G", "g-tol", "Tolerance for gradient in optimization (default: 0 -> function of position tolerance)", false, 0, "gradient tolerance", cmd);
    TCLAP::ValueArg<unsigned int> maxEvalArg("e", "max-eval", "Maximum evaluations (default: 0 -> function of number of unknowns)", false, 0, "max evaluations", cmd);

    TCLAP::ValueArg<unsigned int> nbThreadsArg("T", "nb-threads", "Number of threads to run on (default: all cores)", false, itk::MultiThreader::GetGlobalDefaultNumberOfThreads(), "number of threads", cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef anima::MCMEstimatorImageFilter <float, double> FilterType;
    typedef FilterType::InputImageType InputImageType;
    typedef FilterType::MaskImageType MaskImageType;
    typedef FilterType::CompartmentType CompartmentType;
    typedef FilterType::Pointer FilterPointer;

    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);

    FilterPointer filter = FilterType::New();

    std::cout << "Loading input DWI image..." << std::endl;

    anima::setMultipleImageFilterInputsFromFileName<InputImageType,FilterType>(dwiArg.getValue(), filter);

    // Load gradient table and b-value list
    std::cout << "Importing gradient table and b-values..." << std::endl;

    typedef anima::GradientFileReader < vnl_vector_fixed<double,3>, double > GFReaderType;
    GFReaderType gfReader;
    gfReader.SetGradientFileName(gradsArg.getValue());
    gfReader.SetBValueBaseString(bvalsArg.getValue());
    gfReader.SetGradientIndependentNormalization(bvalueScaleArg.isSet());
    gfReader.Update();

    GFReaderType::GradientVectorType directions = gfReader.GetGradients();
    for(unsigned int i = 0;i < directions.size();++i)
        filter->AddGradientDirection(i,directions[i]);

    GFReaderType::BValueVectorType mb = gfReader.GetBValues();
    filter->SetBValuesList(mb);

    if (computationMaskArg.getValue() != "")
        filter->SetComputationMask(anima::readImage<MaskImageType>(computationMaskArg.getValue()));

    if (inMoseArg.getValue() != "")
        filter->SetMoseVolume(anima::readImage<FilterType::MoseImageType>(inMoseArg.getValue()));

    filter->SetB0Threshold(b0thrArg.getValue());
    
    filter->SetVNLDerivativeComputation(vnlDerivArg.isSet());
    filter->SetAbsoluteInitialDiffusivities(diffInitAbsoluteArg.isSet());

    filter->SetModelWithFreeWaterComponent(freeWaterCompartmentArg.isSet());
    filter->SetModelWithStationaryWaterComponent(stationaryWaterCompartmentArg.isSet());
    filter->SetModelWithRestrictedWaterComponent(restrictedWaterCompartmentArg.isSet());

    switch (compartmentTypeArg.getValue())
    {
        case 1:
            filter->SetCompartmentType(anima::Stick);
            break;

        case 2:
            filter->SetCompartmentType(anima::Zeppelin);
            break;

        case 3:
            filter->SetCompartmentType(anima::Tensor);
            break;
            
        case 4:
            filter->SetCompartmentType(anima::NODDI);
            break;

        default:
            std::cerr << "Unsupported compartment type" << std::endl;
            return EXIT_FAILURE;
    }

    filter->SetNumberOfCompartments(nbFasciclesArg.getValue());
    filter->SetFindOptimalNumberOfCompartments(aicSelectNbCompartmentsArg.isSet());

    filter->SetOptimizer(optimizerArg.getValue());
    filter->SetNumberOfRandomRestarts(nbRandomRestartsArg.getValue());
    filter->SetAbsoluteCostChange(absCostChangeArg.getValue());

    filter->SetNoiseType(FilterType::Gaussian);
    filter->SetMLEstimationStrategy((FilterType::MaximumLikelihoodEstimationMode)mlModeArg.getValue());

    filter->SetXTolerance(xTolArg.getValue());
    filter->SetGTolerance(gTolArg.getValue());
    filter->SetMaxEval(maxEvalArg.getValue());

    filter->SetUseConcentrationBoundsFromDTI(false);
    filter->SetUseFixedWeights(fixWeightsArg.isSet() && (mlModeArg.getValue() != 2));

    filter->SetFreeWaterProportionFixedValue(freeWaterWeightArg.getValue());
    filter->SetStationaryWaterProportionFixedValue(stationaryWaterWeightArg.getValue());
    filter->SetRestrictedWaterProportionFixedValue(restrictedWaterWeightArg.getValue());

    filter->SetUseConstrainedDiffusivity(fixDiffArg.isSet());
    filter->SetUseConstrainedFreeWaterDiffusivity(!optFWDiffArg.isSet());
    filter->SetUseConstrainedIRWDiffusivity(!optIRWDiffArg.isSet());

    if (!fixDiffArg.isSet())
        filter->SetUseCommonDiffusivities(commonDiffusivitiesArg.isSet());
    else
        filter->SetUseCommonDiffusivities(false);

    filter->SetNumberOfThreads(nbThreadsArg.getValue());
    filter->AddObserver(itk::ProgressEvent(), callback);

    itk::TimeProbe tmpTimer;
    tmpTimer.Start();

    try
    {
        filter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }

    tmpTimer.Stop();

    std::cout << "\nEstimation done in " << tmpTimer.GetTotal() << " s" << std::endl;
    std::cout << "Writing MCM to: " << outArg.getValue() << std::endl;

    try
    {
        filter->WriteMCMOutput(outArg.getValue());
    }
    catch (std::exception &e)
    {
        std::cerr << "Error writing output " << e.what() << std::endl;
    }

    if (aicArg.getValue() != "")
    {
        std::cout << "Writing AICu image to: " << aicArg.getValue() << std::endl;
        anima::writeImage(aicArg.getValue(),filter->GetAICcVolume());
    }

    if (outB0Arg.getValue() != "")
    {
        std::cout << "Writing B0 image to: " << outB0Arg.getValue() << std::endl;
        anima::writeImage(outB0Arg.getValue(),filter->GetB0Volume());
    }

    if (outSigmaArg.getValue() != "")
    {
        std::cout << "Writing noise sigma square image to: " << outSigmaArg.getValue() << std::endl;
        anima::writeImage(outSigmaArg.getValue(),filter->GetSigmaSquareVolume());
    }

    if (outMoseArg.getValue() != "")
    {
        std::cout << "Writing model selection image to: " << outMoseArg.getValue() << std::endl;
        anima::writeImage(outMoseArg.getValue(),filter->GetMoseVolume());
    }

    return EXIT_SUCCESS;
}
