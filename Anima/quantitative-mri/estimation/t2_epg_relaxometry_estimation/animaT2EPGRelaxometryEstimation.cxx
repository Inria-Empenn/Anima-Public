#include <iostream>
#include <tclap/CmdLine.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <animaT2EPGRelaxometryEstimationImageFilter.h>
#include <animaReadWriteFunctions.h>
#include <itkTimeProbe.h>

//Update progression of the process
void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
    std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
}

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
	
    TCLAP::ValueArg<std::string> t2Arg("i","t2","List of T2 relaxometry images",true,"","T2 relaxometry images",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask",false,"","computation mask",cmd);

    TCLAP::ValueArg<std::string> t1MapArg("","t1","T1 map",false,"","T1 map",cmd);

    TCLAP::ValueArg<std::string> resT2Arg("o","out-t2","Result T2 image",true,"","result T2 image",cmd);
    TCLAP::ValueArg<std::string> resM0Arg("O","out-m0","Result M0 image",false,"","result M0 image",cmd);
    TCLAP::ValueArg<std::string> resB1Arg("","out-b1","B1 map",false,"","B1 map",cmd);

    TCLAP::SwitchArg nonUniformPulsesArg("N","non-uniform","Use a non uniform pulse profile (default: no)",cmd);
    TCLAP::ValueArg<std::string> excitationProfileArg("E","excitation-profile","Excitation profile text file",false,"","excitation profile file",cmd);
    TCLAP::ValueArg<std::string> pulseProfileArg("p","pulse-profile","Pulse profile text file",false,"","pulse profile file",cmd);
    TCLAP::ValueArg<double> refSliceThicknesshArg("s","ref-thickness","Reference slice thickness used for profile simulation (in mm, default: 3)",false,3,"reference thickness",cmd);
    TCLAP::ValueArg<double> pulseWidthFactorArg("w","pulse-width-factor","Width factor applied at acquisition (default: 1.5)",false,1.5,"pulse width factor",cmd);

    TCLAP::ValueArg<double> trArg("","tr","Repetition time for T2 relaxometry (default: 5000)",false,5000,"repetition time",cmd);
    TCLAP::ValueArg<double> upperBoundT2Arg("u","upper-bound-t2","T2 value upper bound (default: 5000)",false,5000,"T2 value upper bound",cmd);

    TCLAP::ValueArg<double> echoSpacingArg("e","echo-spacing","Spacing between two successive echoes (default: 10)",false,10,"Spacing between echoes",cmd);
    TCLAP::ValueArg<double> excitationT2FlipAngleArg("","t2-ex-flip","Excitation flip angle for T2 (in degrees, default: 90)",false,90,"T2 excitation flip angle",cmd);
    TCLAP::ValueArg<double> t2FlipAngleArg("","t2-flip","All flip angles for T2 (in degrees, default: 180)",false,180,"T2 flip angle",cmd);
    TCLAP::ValueArg<double> backgroundSignalThresholdArg("t","signal-thr","Background signal threshold (default: 10)",false,10,"Background signal threshold",cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("T","numberofthreads","Number of threads to run on (default : all cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
	
    TCLAP::ValueArg<unsigned int> numOptimizerIterArg("","opt-iter","Maximal number of optimizer iterations (default: 2000)",false,2000,"Maximal number of optimizer iterations",cmd);
    TCLAP::ValueArg<double> optimizerStopConditionArg("","opt-stop","Optimizer stopping threshold (default: 1.0e-4)",false,1.0e-4,"Optimizer stopping threshold",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);

    typedef itk::Image <double,3> InputImageType;
    typedef InputImageType OutputImageType;
    typedef anima::T2EPGRelaxometryEstimationImageFilter <InputImageType,OutputImageType> FilterType;
    
    FilterType::Pointer mainFilter = FilterType::New();
	
    // Read and set T2 relaxometry images
    unsigned int numInputs = anima::setMultipleImageFilterInputsFromFileName<InputImageType,FilterType>(t2Arg.getValue(),mainFilter);
   
    mainFilter->SetEchoSpacing(echoSpacingArg.getValue());
    mainFilter->SetT2FlipAngles(t2FlipAngleArg.getValue() * M_PI / 180.0,numInputs);
    mainFilter->SetT2ExcitationFlipAngle(excitationT2FlipAngleArg.getValue() * M_PI / 180.0);

    mainFilter->SetT2UpperBound(upperBoundT2Arg.getValue());
    mainFilter->SetTRValue(trArg.getValue());
    
    mainFilter->SetMaximumOptimizerIterations(numOptimizerIterArg.getValue());
    mainFilter->SetOptimizerStopCondition(optimizerStopConditionArg.getValue());

    mainFilter->SetUniformPulses(!nonUniformPulsesArg.isSet());
    if (nonUniformPulsesArg.isSet())
    {
        mainFilter->SetReferenceSliceThickness(refSliceThicknesshArg.getValue());
        mainFilter->SetPulseWidthFactor(nonUniformPulsesArg.getValue());
        if (pulseProfileArg.getValue() == "")
        {
            std::cerr << "Error: pulse profile needed when using non uniform pulse profiles" << std::endl;
            return EXIT_FAILURE;
        }

        std::vector < std::pair <double, double> > pulseProfile;
        std::ifstream inputPulse(pulseProfileArg.getValue());
        while (!inputPulse.eof())
        {
            char tmpStr[2048];
            inputPulse.getline(tmpStr,2048);

            if (strcmp(tmpStr,"") == 0)
                continue;

            std::stringstream tmpInput;
            tmpInput << tmpStr;

            double xVal, yVal;
            tmpInput >> xVal >> yVal;

            pulseProfile.push_back(std::make_pair(xVal, yVal));
        }

        mainFilter->SetPulseProfile(pulseProfile);

        if (excitationProfileArg.getValue() == "")
        {
            std::cerr << "Error: excitation profile needed when using non uniform pulse profiles" << std::endl;
            return EXIT_FAILURE;
        }

        std::vector < std::pair <double, double> > excitationProfile;
        std::ifstream inputExcitation(excitationProfileArg.getValue());
        while (!inputExcitation.eof())
        {
            char tmpStr[2048];
            inputExcitation.getline(tmpStr,2048);

            if (strcmp(tmpStr,"") == 0)
                continue;

            std::stringstream tmpInput;
            tmpInput << tmpStr;

            double xVal, yVal;
            tmpInput >> xVal >> yVal;

            excitationProfile.push_back(std::make_pair(xVal, yVal));
        }

        inputExcitation.close();
        mainFilter->SetExcitationProfile(excitationProfile);
    }

    if (t1MapArg.getValue() != "")
        mainFilter->SetT1Map(anima::readImage <InputImageType> (t1MapArg.getValue()));
    
    if (maskArg.getValue() != "")
        mainFilter->SetComputationMask(anima::readImage < itk::Image <unsigned char, 3> > (maskArg.getValue()));

    mainFilter->SetAverageSignalThreshold(backgroundSignalThresholdArg.getValue());
    mainFilter->SetNumberOfWorkUnits(nbpArg.getValue());
    
    itk::TimeProbe tmpTime;
    tmpTime.Start();
    
    mainFilter->AddObserver(itk::ProgressEvent(), callback);
    mainFilter->Update();
    
    tmpTime.Stop();
    
    std::cout << "\nTotal computation time: " << tmpTime.GetTotal() << std::endl;

    anima::writeImage <OutputImageType> (resT2Arg.getValue(),mainFilter->GetOutput(0));

    if (resM0Arg.getValue() != "")
        anima::writeImage <OutputImageType> (resM0Arg.getValue(),mainFilter->GetOutput(1));

    if (resB1Arg.getValue() != "")
        anima::writeImage <OutputImageType> (resB1Arg.getValue(),mainFilter->GetOutput(2));

    return EXIT_SUCCESS;
}
