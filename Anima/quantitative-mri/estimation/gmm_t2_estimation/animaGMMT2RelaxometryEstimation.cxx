#include <iostream>
#include <tclap/CmdLine.h>

#include <animaGMMT2RelaxometryEstimationImageFilter.h>
#include <animaReadWriteFunctions.h>
#include <itkTimeProbe.h>

#include <itkCommand.h>

//Update progression of the process
void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
    std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
}

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
	
    TCLAP::ValueArg<std::string> t2Arg("i","t2","T2 relaxometry images (list or 4D volume)",true,"","T2 relaxometry images",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask",false,"","computation mask",cmd);

    TCLAP::ValueArg<std::string> t1MapArg("","t1","T1 map",false,"","T1 map",cmd);
    TCLAP::ValueArg<std::string> gaussMeansArg("g","gauss-mean","Text file with Gaussian means",false,"","Gaussian means",cmd);
    TCLAP::ValueArg<std::string> gaussVarsArg("G","gauss-var","Text file with Gaussian variances",false,"","Gaussian variances",cmd);

    TCLAP::ValueArg<std::string> resWeightsArg("O","out-weights","Result weights image",false,"","result weights image",cmd);
    TCLAP::ValueArg<std::string> resM0Arg("","out-m0","Result M0 image",false,"","result M0 image",cmd);
    TCLAP::ValueArg<std::string> resMWFArg("o","out-mwf","Result MWF image",true,"","result MWF image",cmd);
    TCLAP::ValueArg<std::string> resB1Arg("","out-b1","Result B1 image",false,"","result B1 image",cmd);
    TCLAP::ValueArg<std::string> resSigmaSqArg("","out-sig","Result sigma square image",false,"","result sigma square image",cmd);

    TCLAP::SwitchArg nonUniformPulsesArg("N","non-uniform","Use a non uniform pulse profile (default: no)",cmd);
    TCLAP::ValueArg<std::string> excitationProfileArg("E","excitation-profile","Excitation profile text file",false,"","excitation profile file",cmd);
    TCLAP::ValueArg<std::string> pulseProfileArg("p","pulse-profile","Pulse profile text file",false,"","pulse profile file",cmd);
    TCLAP::ValueArg<double> refSliceThicknesshArg("s","ref-thickness","Reference slice thickness used for profile simulation (in mm, default: 3)",false,3,"reference thickness",cmd);
    TCLAP::ValueArg<double> pulseWidthFactorArg("w","pulse-width-factor","Width factor applied at acquisition (default: 1.5)",false,1.5,"pulse width factor",cmd);

    TCLAP::ValueArg<double> echoSpacingArg("e","echo-spacing","Spacing between two successive echoes (default: 10)",false,10,"Spacing between echoes",cmd);
    TCLAP::ValueArg<double> excitationT2FlipAngleArg("","t2-ex-flip","Excitation flip angle for T2 (in degrees, default: 90)",false,90,"T2 excitation flip angle",cmd);
    TCLAP::ValueArg<double> t2FlipAngleArg("","t2-flip","All flip angles for T2 (in degrees, default: 180)",false,180,"T2 flip angle",cmd);
    TCLAP::ValueArg<double> backgroundSignalThresholdArg("t","signal-thr","Background signal threshold (default: 10)",false,10,"Background signal threshold",cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("T","numberofthreads","Number of threads to run on (default : all cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
	
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
    typedef anima::GMMT2RelaxometryEstimationImageFilter <double> FilterType;
    typedef FilterType::OutputImageType OutputImageType;
    typedef FilterType::VectorOutputImageType VectorOutputImageType;
    typedef itk::ImageFileReader <InputImageType> InputImageReaderType;

    FilterType::Pointer mainFilter = FilterType::New();
	
    unsigned int numInputs = anima::setMultipleImageFilterInputsFromFileName<InputImageType,FilterType>(t2Arg.getValue(),mainFilter);

    mainFilter->SetEchoSpacing(echoSpacingArg.getValue());
    mainFilter->SetT2FlipAngles(t2FlipAngleArg.getValue() * M_PI / 180.0,numInputs);
    mainFilter->SetT2ExcitationFlipAngle(excitationT2FlipAngleArg.getValue() * M_PI / 180.0);

    mainFilter->SetUniformPulses(!nonUniformPulsesArg.isSet());
    if (nonUniformPulsesArg.isSet())
    {
        mainFilter->SetReferenceSliceThickness(refSliceThicknesshArg.getValue());
        mainFilter->SetPulseWidthFactor(pulseWidthFactorArg.getValue());
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

    if (gaussMeansArg.getValue() != "")
        mainFilter->SetGaussianMeans(gaussMeansArg.getValue());

    if (gaussVarsArg.getValue() != "")
        mainFilter->SetGaussianVariances(gaussVarsArg.getValue());

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
    mainFilter->SetNumberOfWorkUnits(nbpArg.getValue());

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
        return EXIT_FAILURE;
    }

    tmpTime.Stop();
    
    std::cout << "\nRegular estimation computation time: " << tmpTime.GetTotal() << std::endl;

    anima::writeImage<OutputImageType> (resMWFArg.getValue(),mainFilter->GetMWFOutputImage());

    if (resWeightsArg.getValue() != "")
        anima::writeImage<VectorOutputImageType> (resWeightsArg.getValue(),mainFilter->GetWeightsImage());

    if (resM0Arg.getValue() != "")
        anima::writeImage<OutputImageType> (resM0Arg.getValue(),mainFilter->GetM0OutputImage());

    if (resB1Arg.getValue() != "")
        anima::writeImage<OutputImageType> (resB1Arg.getValue(),mainFilter->GetB1OutputImage());

    if (resSigmaSqArg.getValue() != "")
        anima::writeImage<OutputImageType> (resSigmaSqArg.getValue(),mainFilter->GetSigmaSquareOutputImage());

    return EXIT_SUCCESS;
}
