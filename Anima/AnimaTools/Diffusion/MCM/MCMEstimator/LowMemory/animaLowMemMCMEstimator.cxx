#include <iostream>
#include <tclap/CmdLine.h>

#include "animaLowMemMCMEstimatorBridge.h"

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    // Required arguments
    TCLAP::ValueArg<std::string> dwiArg("i", "dwi-list", "DWI volume list", true, "", "DWI images", cmd);
    TCLAP::ValueArg<std::string> gradsArg("g", "grads", "Gradient table", true, "", "gradients", cmd);
    TCLAP::ValueArg<std::string> bvalsArg("b", "bvals", "B-value list", true, "", "b-values", cmd);
    TCLAP::ValueArg<std::string> inMoseArg("", "in-mose", "Input model selection map (overrides model selection step)", false, "", "model selection input image", cmd);
    TCLAP::SwitchArg bvalueScaleArg("B","b-no-scale","Do not scale b-values according to gradient norm",cmd);
    TCLAP::ValueArg<double> smallDeltaArg("", "small-delta", "Diffusion small delta (in seconds)", false, anima::DiffusionSmallDelta, "small delta", cmd);
    TCLAP::ValueArg<double> bigDeltaArg("", "big-delta", "Diffusion big delta (in seconds)", false, anima::DiffusionBigDelta, "big delta", cmd);

    // Outputs
    TCLAP::ValueArg<std::string> outArg("o", "out-mcm", "MCM output prefix", true, "", "MCM output prefix", cmd);
    TCLAP::ValueArg<std::string> aicArg("a", "out-aic", "Output estimated AICu prefix", false, "", "AICu output prefix", cmd);
    TCLAP::ValueArg<std::string> outB0Arg("", "out-b0", "Output estimated B0 prefix", false, "", "B0 output prefix", cmd);
    TCLAP::ValueArg<std::string> outSigmaArg("", "out-sig", "Output estimated noise sigma square prefix", false, "", "noise sigma square output prefix", cmd);
    TCLAP::ValueArg<std::string> outMoseArg("", "out-mose", "Output model selection map prefix", false, "", "model selection output prefix", cmd);

    // Optional arguments
    TCLAP::ValueArg<std::string> computationMaskArg("m", "mask", "Computation mask", false, "", "computation mask", cmd);
    TCLAP::ValueArg<double> b0thrArg("", "b0-thr", "Background threshold on B0 value (default: 10)", false, 10.0, "B0 theshold", cmd);

    TCLAP::ValueArg<unsigned int> nbFasciclesArg("n", "nb-fascicles", "Number of computed fascicles (default: 2)", false, 2, "number of fascicles", cmd);
    TCLAP::ValueArg<unsigned int> compartmentTypeArg("c", "comp-type", "Compartment type for fascicles: 1: stick, 2: zeppelin, 3: tensor, 4: NODDI, 5: DDI (default: 3)", false, 3, "fascicles type", cmd);
    TCLAP::SwitchArg aicSelectNbCompartmentsArg("M", "opt-nb-comp", "Activate AICC-based number of compartments selection", cmd, false);

    TCLAP::SwitchArg freeWaterCompartmentArg("F", "free-water", "Model with free water", cmd, false);
    TCLAP::SwitchArg stationaryWaterCompartmentArg("W", "stationary-water", "Model with stationary water", cmd, false);
    TCLAP::SwitchArg restrictedWaterCompartmentArg("R", "restricted-water", "Model with restricted water", cmd, false);
    TCLAP::SwitchArg staniszCompartmentArg("Z", "stanisz", "Model with stanisz isotropic compartment", cmd, false);

    TCLAP::SwitchArg optFWDiffArg("", "opt-free-water-diff", "Optimize free water diffusivity value", cmd, false);
    TCLAP::SwitchArg optIRWDiffArg("", "opt-ir-water-diff", "Optimize isotropic restricted water diffusivity value", cmd, false);
    TCLAP::SwitchArg optStaniszRadiusArg("", "opt-stanisz-radius", "Optimize isotropic Stanisz radius value", cmd, false);
    TCLAP::SwitchArg optStaniszDiffArg("", "opt-stanisz-diff", "Optimize isotropic Stanisz diffusivity value", cmd, false);

    TCLAP::SwitchArg fixDiffArg("", "fix-diff", "Fix diffusivity value", cmd, false);
    TCLAP::SwitchArg fixKappaArg("", "fix-kappa", "Fix orientation concentration values", cmd, false);
    TCLAP::SwitchArg fixEAFArg("", "fix-eaf", "Fix extra axonal fraction values", cmd, false);

    TCLAP::SwitchArg commonDiffusivitiesArg("", "common-diffusivities", "Share diffusivity values among compartments", cmd, false);
    TCLAP::SwitchArg commonKappaArg("", "common-kappa", "Share orientation concentration values among compartments", cmd, false);
    TCLAP::SwitchArg commonEAFArg("", "common-eaf", "Share extra axonal fraction values among compartments", cmd, false);

    //Initial values for diffusivities
    TCLAP::ValueArg<double> initAxialDiffArg("", "init-axial-diff", "Initial axial diffusivity (default: 1.71e-3)", false, 1.71e-3, "initial axial diffusivity", cmd);
    TCLAP::ValueArg<double> initRadialDiff1Arg("", "init-radial-diff1", "Initial first radial diffusivity (default: 1.9e-4)", false, 1.9e-4, "initial first radial diffusivity", cmd);
    TCLAP::ValueArg<double> initRadialDiff2Arg("", "init-radial-diff2", "Initial second radial diffusivity (default: 1.5e-4)", false, 1.5e-4, "initial second radial diffusivity", cmd);
    TCLAP::ValueArg<double> initIRWDiffArg("", "init-irw-diff", "Initial isotropic restricted diffusivity (default: 7.5e-4)", false, 7.5e-4, "initial IRW diffusivity", cmd);
    TCLAP::ValueArg<double> initStaniszDiffArg("", "init-stanisz-diff", "Initial Stanisz diffusivity (default: 1.71e-3)", false, 1.71e-3, "initial Stanisz diffusivity", cmd);

    // Optimization parameters
    TCLAP::ValueArg<std::string> optimizerArg("", "optimizer", "Optimizer for estimation: bobyqa (default), ccsaq, bfgs or levenberg", false, "bobyqa", "optimizer", cmd);
    TCLAP::ValueArg<double> absCostChangeArg("", "abs-cost-change", "Cost function change to stop estimation (default: 0.01)", false, 0.01, "cost change threshold", cmd);
    TCLAP::ValueArg <unsigned int> mlModeArg("", "ml-mode", "ML estimation strategy: marginal likelihood (0), profile likelihood (1, default), Variable projection (2)", false, 1, "ML mode", cmd);
    TCLAP::ValueArg<double> xTolArg("x", "x-tol", "Tolerance for position in optimization (default: 0 -> 1.0e-4 or 1.0e-7 for bobyqa)", false, 0, "position tolerance", cmd);
    TCLAP::ValueArg<double> fTolArg("", "f-tol", "Tolerance for relative cost in optimization (default: 0 -> function of position tolerance)", false, 0, "cost relative tolerance", cmd);
    TCLAP::ValueArg<unsigned int> maxEvalArg("e", "max-eval", "Maximum evaluations (default: 0 -> function of number of unknowns)", false, 0, "max evaluations", cmd);

    TCLAP::ValueArg<unsigned int> noiseTypeArg("", "noise-type", "Noise type for optimization: 0: Gaussian (default) or 1: NCC", false, 0, "noise type", cmd);
    TCLAP::ValueArg<unsigned int> nbCoilsArg("", "nb-coils", "Number of coils, for NCC estimation (default: 1)", false, 1, "number of coils", cmd);

    TCLAP::ValueArg<unsigned int> nbThreadsArg("T", "nb-threads", "Number of threads to run on (default: all cores)", false, itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(), "number of threads", cmd);

    TCLAP::ValueArg<unsigned int> splitsArg("s","split","Split image for low memory (default: 2)",false,2,"Number of splits",cmd);
    TCLAP::ValueArg<int> specSplitArg("S","splittoprocess","Specific split to process (use to run on cluster (default: -1 = all)",false,-1,"Split to process",cmd);
    TCLAP::SwitchArg genOutputDescroArg("G","generateouputdescription","Generate ouptut description data",cmd,false);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }
    	
    typedef anima::LowMemMCMEstimatorBridge EstimatorBridgeType;
    
    EstimatorBridgeType *mainFilter = new EstimatorBridgeType;

    mainFilter->SetComputationMask(computationMaskArg.getValue());
    mainFilter->SetNumberOfThreads(nbThreadsArg.getValue());

    mainFilter->SetDWIFileNames(dwiArg.getValue());
    if (inMoseArg.getValue() != "")
        mainFilter->SetInputMoseName(inMoseArg.getValue());

    mainFilter->SetOutputName(outArg.getValue());
	mainFilter->SetNbSplits(splitsArg.getValue());

    mainFilter->SetGradients(gradsArg.getValue());
    mainFilter->SetBValues(bvalsArg.getValue());
    mainFilter->SetBValueScale(bvalueScaleArg.getValue());

    mainFilter->SetOutputAICName(aicArg.getValue());
    mainFilter->SetOutputB0Name(outB0Arg.getValue());
    mainFilter->SetOutputSigmaName(outSigmaArg.getValue());
    mainFilter->SetOutputMoseName(outMoseArg.getValue());

    mainFilter->SetB0Threshold(b0thrArg.getValue());
    mainFilter->SetNumberOfFascicles(nbFasciclesArg.getValue());
    mainFilter->SetFindOptimalNumberOfCompartments(aicSelectNbCompartmentsArg.isSet());
    mainFilter->SetCompartmentType(compartmentTypeArg.getValue());

    mainFilter->SetFreeWaterCompartment(freeWaterCompartmentArg.isSet());
    mainFilter->SetStationaryWaterCompartment(stationaryWaterCompartmentArg.isSet());
    mainFilter->SetRestrictedWaterCompartment(restrictedWaterCompartmentArg.isSet());
    mainFilter->SetStaniszCompartment(staniszCompartmentArg.isSet());

    mainFilter->SetAxialDiffusivityValue(initAxialDiffArg.getValue());
    mainFilter->SetRadialDiffusivity1Value(initRadialDiff1Arg.getValue());
    mainFilter->SetRadialDiffusivity2Value(initRadialDiff2Arg.getValue());
    mainFilter->SetIRWDiffusivityValue(initIRWDiffArg.getValue());
    mainFilter->SetStaniszDiffusivityValue(initStaniszDiffArg.getValue());

    mainFilter->SetOptimizeFreeWaterDiffusivity(optFWDiffArg.isSet());
    mainFilter->SetOptimizeIRWDiffusivity(optIRWDiffArg.isSet());
    mainFilter->SetOptimizeStaniszDiffusivity(optStaniszDiffArg.isSet());
    mainFilter->SetOptimizeStaniszRadius(optStaniszRadiusArg.isSet());

    mainFilter->SetFixDiffusivity(fixDiffArg.isSet());
    mainFilter->SetFixKappa(fixKappaArg.isSet());
    mainFilter->SetFixEAF(fixEAFArg.isSet());

    mainFilter->SetCommonDiffusivities(commonDiffusivitiesArg.isSet());
    mainFilter->SetCommonKappa(commonKappaArg.isSet());
    mainFilter->SetCommonEAF(commonEAFArg.isSet());

    mainFilter->SetOptimizerType(optimizerArg.getValue());
    mainFilter->SetAbsCostChange(absCostChangeArg.getValue());

    mainFilter->SetXTolerance(xTolArg.getValue());
    mainFilter->SetFTolerance(fTolArg.getValue());
    mainFilter->SetMaxEval(maxEvalArg.getValue());

    mainFilter->SetNoiseType((EstimatorBridgeType::SignalNoiseType)noiseTypeArg.getValue());
    mainFilter->SetNumberOfCoils(nbCoilsArg.getValue());
    mainFilter->SetMLEstimationStrategy((EstimatorBridgeType::MaximumLikelihoodEstimationMode)mlModeArg.getValue());

    mainFilter->SetSmallDelta(smallDeltaArg.getValue());
    mainFilter->SetBigDelta(bigDeltaArg.getValue());

    try
    {
        mainFilter->Update(specSplitArg.getValue(),genOutputDescroArg.getValue());
    }
    catch(itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }
    
	delete mainFilter;
	
    return EXIT_SUCCESS;
}
