#include <iostream>
#include <tclap/CmdLine.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <animaMultiT2RelaxometryEstimationImageFilter.h>
#include <animaReadWriteFunctions.h>
#include <itkTimeProbe.h>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <itkCommand.h>

//Update progression of the process
void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
    std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
}

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);
	
    TCLAP::ValueArg<std::string> t2Arg("l","t2","List of T2 relaxometry images",true,"","T2 relaxometry images",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask",false,"","computation mask",cmd);

    TCLAP::ValueArg<std::string> t1MapArg("","t1","T1 map",false,"","T1 map",cmd);

    TCLAP::ValueArg<std::string> resT2Arg("O","out-t2","Result T2 image",false,"","result T2 image",cmd);
    TCLAP::ValueArg<std::string> resM0Arg("","out-m0","Result M0 image",false,"","result M0 image",cmd);
    TCLAP::ValueArg<std::string> resMWFArg("o","out-mwf","Result MWF image",true,"","result MWF image",cmd);
    TCLAP::ValueArg<std::string> resB1Arg("","out-b1","Result B1 image",false,"","result B1 image",cmd);
    TCLAP::ValueArg<std::string> resCostArg("c","out-cost","Result cost image",false,"","result cost image",cmd);

    TCLAP::ValueArg<double> echoSpacingArg("e","echo-spacing","Spacing between two successive echoes (default: 10)",false,10,"Spacing between echoes",cmd);
    TCLAP::ValueArg<double> excitationT2FlipAngleArg("","t2-ex-flip","Excitation flip angle for T2 (in degrees, default: 90)",false,90,"T2 excitation flip angle",cmd);
    TCLAP::ValueArg<double> t2FlipAngleArg("","t2-flip","All flip angles for T2 (in degrees, default: 180)",false,180,"T2 flip angle",cmd);
    TCLAP::ValueArg<double> backgroundSignalThresholdArg("t","signal-thr","Background signal threshold (default: 10)",false,10,"Background signal threshold",cmd);
    TCLAP::ValueArg<double> regulIntensityArg("r","reg-int","Tikhonov regularization intensity (default: 1.08)",false,1.08,"regularization intensity",cmd);
    TCLAP::ValueArg<double> nlRegulIntensityArg("","nl-reg-int","Tikhonov regularization intensity in NL mode (default: 1.5)",false,1.5,"regularization intensity",cmd);

    TCLAP::ValueArg<double> myelinThrArg("","mye-thr","T2 myelin threshold for MWF computation (default: 50)",false,50,"T2 myelin threshold",cmd);
    TCLAP::ValueArg<unsigned int> numT2CompartmentsArg("n","num-t2","Number of T2 compartments (default: 40)",false,40,"Number of T2 compartments",cmd);
    TCLAP::ValueArg<double> t2LowerBoundArg("","low-t2","Lower T2 value (default: 15)",false,15,"T2 lower value",cmd);
    TCLAP::ValueArg<double> t2UpperBoundArg("","up-t2","Upper T2 value (default: 2000)",false,2000,"T2 upper value",cmd);

    TCLAP::SwitchArg b1OnExcAngleArg("B", "b1-exc", "B1 is also applied to excitation angle",cmd,false);
    TCLAP::ValueArg<unsigned int> b1NumOptimizerIterArg("","b1-opt-iter","Maximal number of optimizer iterations (default: 200)",false,200,"Maximal number of optimizer iterations",cmd);
    TCLAP::ValueArg<double> b1OptimizerStopConditionArg("","b1-opt-stop","Optimizer stopping threshold (default: 1.0e-4)",false,1.0e-4,"B1 Optimizer stopping threshold",cmd);
    TCLAP::ValueArg<double> b1OptimizerInitialStepArg("i","b1-opt-init","Optimizer initial step (default: 0.1)",false,0.1,"B1 Optimizer initial step",cmd);
    TCLAP::ValueArg<double> b1ToleranceArg("","b1-tol","B1 tolerance threshold of convergence (default: 1.0e-4)",false,1.0e-4,"B1 tolerance for optimization convergence",cmd);

    //NL params
    TCLAP::SwitchArg nlEstimationArg("N","nl-est","NL estimation of MWF",cmd, false);
    TCLAP::ValueArg<double> weightThrArg("w","weightThr","Weight threshold: patches around have to be similar enough -> default: 0.0",false,0.0,"Weight threshold",cmd);
    TCLAP::ValueArg<double> betaArg("b","beta","Beta parameter for local noise estimation -> default: 1",false,1,"Beta for local noise estimation",cmd);
    TCLAP::ValueArg<double> meanMinArg("","meanMin","Minimun mean threshold (default: 0.95)",false,0.95,"Minimun mean threshold",cmd);
    TCLAP::ValueArg<double> varMinArg("v","varMin","Minimun variance threshold -> default: 0.5",false,0.5,"Minimun variance threshold",cmd);
    TCLAP::ValueArg<unsigned int> patchHSArg("S","patchHalfSize","Patch half size in each direction -> default: 1",false,1,"patch half size",cmd);
    TCLAP::ValueArg<unsigned int> patchSSArg("s","patchStepSize","Patch step size for searching -> default: 1",false,1,"Patch search step size",cmd);
    TCLAP::ValueArg<unsigned int> patchNeighArg("","patchNeighborhood","Patch half neighborhood size -> default: 5",false,5,"Patch search neighborhood size",cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("T","numberofthreads","Number of threads to run on (default : all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
	
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }
    
    typedef itk::Image <double,3> InputImageType;
    typedef InputImageType OutputImageType;
    typedef anima::MultiT2RelaxometryEstimationImageFilter <double> FilterType;
    typedef FilterType::OutputImageType OutputImageType;
    typedef FilterType::VectorOutputImageType VectorOutputImageType;
    typedef itk::ImageFileReader <InputImageType> InputImageReaderType;

    FilterType::Pointer mainFilter = FilterType::New();
	
    unsigned int numInputs = anima::setMultipleImageFilterInputsFromFileName<InputImageType,FilterType>(t2Arg.getValue(),mainFilter);

    mainFilter->SetEchoSpacing(echoSpacingArg.getValue());
    mainFilter->SetT2FlipAngles(t2FlipAngleArg.getValue() * M_PI / 180.0,numInputs);
    mainFilter->SetT2ExcitationFlipAngle(excitationT2FlipAngleArg.getValue() * M_PI / 180.0);
    mainFilter->SetB1OnExcitationAngle(b1OnExcAngleArg.isSet());
    mainFilter->SetB1MaximumOptimizerIterations(b1NumOptimizerIterArg.getValue());
    mainFilter->SetB1OptimizerStopCondition(b1OptimizerStopConditionArg.getValue());
    mainFilter->SetB1OptimizerInitialStep(b1OptimizerInitialStepArg.getValue());
    mainFilter->SetB1Tolerance(b1ToleranceArg.getValue());

    mainFilter->SetLowerT2Bound(t2LowerBoundArg.getValue());
    mainFilter->SetUpperT2Bound(t2UpperBoundArg.getValue());
    mainFilter->SetMyelinThreshold(myelinThrArg.getValue());
    mainFilter->SetNumberOfT2Compartments(numT2CompartmentsArg.getValue());
    mainFilter->SetRegularizationIntensity(regulIntensityArg.getValue());
    mainFilter->SetNLEstimation(false);

    if (t1MapArg.getValue() != "")
    {
        InputImageReaderType::Pointer t1MapRead = InputImageReaderType::New();
        t1MapRead->SetFileName(t1MapArg.getValue().c_str());
        t1MapRead->Update();
        
        mainFilter->SetT1Map(t1MapRead->GetOutput());
    }
    
    if (maskArg.getValue() != "")
    {
        typedef itk::ImageFileReader < itk::Image <unsigned char, 3> > itkMaskReader;
        itkMaskReader::Pointer maskRead = itkMaskReader::New();
        maskRead->SetFileName(maskArg.getValue().c_str());
        maskRead->Update();
        
        mainFilter->SetComputationMask(maskRead->GetOutput());
    }
    
    mainFilter->SetAverageSignalThreshold(backgroundSignalThresholdArg.getValue());
    mainFilter->SetNumberOfThreads(nbpArg.getValue());

    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);
    mainFilter->AddObserver(itk::ProgressEvent(), callback );

    itk::TimeProbe tmpTime;
    tmpTime.Start();
    
    try
    {
        mainFilter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
    }

    tmpTime.Stop();
    
    std::cout << "Regular estimation computation time: " << tmpTime.GetTotal() << std::endl;

    FilterType::OutputImagePointer outMWFImage = mainFilter->GetMWFOutputImage();
    outMWFImage->DisconnectPipeline();
    FilterType::OutputImagePointer outM0Image = mainFilter->GetM0OutputImage();
    outM0Image->DisconnectPipeline();
    FilterType::OutputImagePointer outCostImage = mainFilter->GetCostOutputImage();
    outCostImage->DisconnectPipeline();
    FilterType::OutputImagePointer outB1Image = mainFilter->GetB1OutputImage();
    outB1Image->DisconnectPipeline();
    FilterType::VectorOutputImagePointer outT2Image = mainFilter->GetT2OutputImage();
    outT2Image->DisconnectPipeline();

    if (nlEstimationArg.isSet())
    {
        FilterType::Pointer secondaryFilter = FilterType::New();
        for (unsigned int i = 0;i < mainFilter->GetNumberOfIndexedInputs();++i)
            secondaryFilter->SetInput(i,mainFilter->GetInput(i));

        secondaryFilter->SetEchoSpacing(echoSpacingArg.getValue());
        secondaryFilter->SetT2FlipAngles(t2FlipAngleArg.getValue() * M_PI / 180.0,numInputs);
        secondaryFilter->SetT2ExcitationFlipAngle(excitationT2FlipAngleArg.getValue() * M_PI / 180.0);
        secondaryFilter->SetB1OnExcitationAngle(b1OnExcAngleArg.isSet());
        secondaryFilter->SetB1MaximumOptimizerIterations(b1NumOptimizerIterArg.getValue());
        secondaryFilter->SetB1OptimizerStopCondition(b1OptimizerStopConditionArg.getValue());
        secondaryFilter->SetB1OptimizerInitialStep(b1OptimizerInitialStepArg.getValue());
        secondaryFilter->SetB1Tolerance(b1ToleranceArg.getValue());

        secondaryFilter->SetLowerT2Bound(t2LowerBoundArg.getValue());
        secondaryFilter->SetUpperT2Bound(t2UpperBoundArg.getValue());
        secondaryFilter->SetMyelinThreshold(myelinThrArg.getValue());
        secondaryFilter->SetNumberOfT2Compartments(numT2CompartmentsArg.getValue());
        secondaryFilter->SetRegularizationIntensity(nlRegulIntensityArg.getValue());

        secondaryFilter->SetT1Map(mainFilter->GetT1Map());
        secondaryFilter->SetComputationMask(mainFilter->GetComputationMask());

        secondaryFilter->SetAverageSignalThreshold(backgroundSignalThresholdArg.getValue());
        secondaryFilter->SetNumberOfThreads(nbpArg.getValue());

        itk::CStyleCommand::Pointer secondaryCallback = itk::CStyleCommand::New();
        secondaryCallback->SetCallback(eventCallback);
        secondaryFilter->AddObserver(itk::ProgressEvent(), secondaryCallback);

        secondaryFilter->SetNLEstimation(true);
        secondaryFilter->SetWeightThreshold(weightThrArg.getValue());
        secondaryFilter->SetBetaParameter(betaArg.getValue());
        secondaryFilter->SetMeanMinThreshold(meanMinArg.getValue());
        secondaryFilter->SetVarMinThreshold(varMinArg.getValue());
        secondaryFilter->SetPatchHalfSize(patchHSArg.getValue());
        secondaryFilter->SetSearchStepSize(patchSSArg.getValue());
        secondaryFilter->SetSearchNeighborhood(patchNeighArg.getValue());

        secondaryFilter->SetInitialB1Map(outB1Image);
        secondaryFilter->SetInitialT2Map(outT2Image);
        secondaryFilter->SetInitialM0Map(outM0Image);

        itk::TimeProbe tmpTimeNL;
        tmpTimeNL.Start();

        try
        {
            secondaryFilter->Update();
        }
        catch (itk::ExceptionObject &e)
        {
            std::cerr << e << std::endl;
        }

        tmpTimeNL.Stop();

        std::cout << "NL estimation computation time: " << tmpTimeNL.GetTotal() << std::endl;

        outMWFImage = secondaryFilter->GetMWFOutputImage();
        outMWFImage->DisconnectPipeline();
        outM0Image = secondaryFilter->GetM0OutputImage();
        outM0Image->DisconnectPipeline();
        outCostImage = secondaryFilter->GetCostOutputImage();
        outCostImage->DisconnectPipeline();
        outB1Image = secondaryFilter->GetB1OutputImage();
        outB1Image->DisconnectPipeline();
        outT2Image = secondaryFilter->GetT2OutputImage();
        outT2Image->DisconnectPipeline();
    }

    anima::writeImage<OutputImageType> (resMWFArg.getValue(),outMWFImage);

    if (resT2Arg.getValue() != "")
        anima::writeImage<VectorOutputImageType> (resT2Arg.getValue(),outT2Image);

    if (resM0Arg.getValue() != "")
        anima::writeImage<OutputImageType> (resM0Arg.getValue(),outM0Image);

    if (resB1Arg.getValue() != "")
        anima::writeImage<OutputImageType> (resB1Arg.getValue(),outB1Image);

    if (resCostArg.getValue() != "")
        anima::writeImage<OutputImageType> (resCostArg.getValue(),outCostImage);

    return EXIT_SUCCESS;
}
