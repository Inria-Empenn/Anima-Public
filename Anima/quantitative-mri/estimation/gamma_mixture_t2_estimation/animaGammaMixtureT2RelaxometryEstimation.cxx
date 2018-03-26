#include <iostream>
#include <tclap/CmdLine.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <animaGammaMixtureT2RelaxometryEstimationImageFilter.h>
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
	
    TCLAP::ValueArg<std::string> t2Arg("i","t2","T2 relaxometry images (list or 4D volume)",true,"","T2 relaxometry images",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask",false,"","computation mask",cmd);

    TCLAP::ValueArg<std::string> t1MapArg("","t1","T1 map",false,"","T1 map",cmd);
    TCLAP::SwitchArg constrainArg("C","constrain","Constrain estimation of means to just the central one",cmd);

    TCLAP::ValueArg<std::string> resWeightsArg("O","out-weights","Result weights image",false,"","result weights image",cmd);
    TCLAP::ValueArg<std::string> resM0Arg("","out-m0","Result M0 image",false,"","result M0 image",cmd);
    TCLAP::ValueArg<std::string> resMWFArg("o","out-mwf","Result MWF image",true,"","result MWF image",cmd);
    TCLAP::ValueArg<std::string> resB1Arg("","out-b1","Result B1 image",false,"","result B1 image",cmd);
    TCLAP::ValueArg<std::string> resMeanParamArg("","mean-param","Result mean parameter image",false,"","result mean parameter image",cmd);

    TCLAP::ValueArg<double> echoSpacingArg("e","echo-spacing","Spacing between two successive echoes (default: 10)",false,10,"Spacing between echoes",cmd);
    TCLAP::ValueArg<double> excitationT2FlipAngleArg("","t2-ex-flip","Excitation flip angle for T2 (in degrees, default: 90)",false,90,"T2 excitation flip angle",cmd);
    TCLAP::ValueArg<double> t2FlipAngleArg("","t2-flip","All flip angles for T2 (in degrees, default: 180)",false,180,"T2 flip angle",cmd);
    TCLAP::ValueArg<double> backgroundSignalThresholdArg("t","signal-thr","Background signal threshold (default: 10)",false,10,"Background signal threshold",cmd);

    TCLAP::ValueArg<double> t2IntegrationStepArg("","t2-int","T2 distribution integration step (default: 1)",false,1.0,"T2 integration step",cmd);
    TCLAP::ValueArg<double> t2LowerBoundArg("","low-t2","Lower T2 value (default: 0)",false,0,"T2 lower value",cmd);
    TCLAP::ValueArg<double> t2UpperBoundArg("","up-t2","Upper T2 value (default: 3000)",false,3000,"T2 upper value",cmd);

    TCLAP::SwitchArg b1OnExcAngleArg("B", "b1-exc", "B1 is also applied to excitation angle",cmd,false);
    TCLAP::ValueArg<double> b1ToleranceArg("","b1-tol","B1 tolerance threshold of convergence (default: 1.0e-4)",false,1.0e-4,"B1 tolerance for optimization convergence",cmd);
    TCLAP::ValueArg<double> costToleranceArg("c","cost-tol","Cost tolerance threshold of convergence (default: 1.0e-4)",false,1.0e-4,"Cost tolerance for optimization convergence",cmd);

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
    typedef anima::GammaMixtureT2RelaxometryEstimationImageFilter <double> FilterType;
    typedef FilterType::OutputImageType OutputImageType;
    typedef FilterType::VectorOutputImageType VectorOutputImageType;
    typedef itk::ImageFileReader <InputImageType> InputImageReaderType;

    FilterType::Pointer mainFilter = FilterType::New();
	
    unsigned int numInputs = anima::setMultipleImageFilterInputsFromFileName<InputImageType,FilterType>(t2Arg.getValue(),mainFilter);

    mainFilter->SetEchoSpacing(echoSpacingArg.getValue());
    mainFilter->SetT2FlipAngles(t2FlipAngleArg.getValue() * M_PI / 180.0,numInputs);
    mainFilter->SetT2ExcitationFlipAngle(excitationT2FlipAngleArg.getValue() * M_PI / 180.0);
    mainFilter->SetB1OnExcitationAngle(b1OnExcAngleArg.isSet());
    mainFilter->SetB1Tolerance(b1ToleranceArg.getValue());
    mainFilter->SetCostTolerance(costToleranceArg.getValue());
    mainFilter->SetConstrainedParameters(constrainArg.isSet());

    mainFilter->SetLowerT2Bound(t2LowerBoundArg.getValue());
    mainFilter->SetUpperT2Bound(t2UpperBoundArg.getValue());
    mainFilter->SetT2IntegrationStep(t2IntegrationStepArg.getValue());

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
    
    std::cout << "Estimation computation time: " << tmpTime.GetTotal() << std::endl;

    anima::writeImage<OutputImageType> (resMWFArg.getValue(),mainFilter->GetMWFOutputImage());

    if (resWeightsArg.getValue() != "")
        anima::writeImage<VectorOutputImageType> (resWeightsArg.getValue(),mainFilter->GetWeightsImage());
	
    if (resMeanParamArg.getValue() != "")
        anima::writeImage<VectorOutputImageType> (resMeanParamArg.getValue(),mainFilter->GetMeanParamImage());

    if (resM0Arg.getValue() != "")
        anima::writeImage<OutputImageType> (resM0Arg.getValue(),mainFilter->GetM0OutputImage());

    if (resB1Arg.getValue() != "")
        anima::writeImage<OutputImageType> (resB1Arg.getValue(),mainFilter->GetB1OutputImage());

    return EXIT_SUCCESS;
}
