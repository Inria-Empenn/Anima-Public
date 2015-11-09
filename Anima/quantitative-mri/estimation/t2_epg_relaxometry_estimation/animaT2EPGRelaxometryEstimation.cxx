#include <iostream>
#include <tclap/CmdLine.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <animaT2EPGRelaxometryEstimationImageFilter.h>
#include <animaReadWriteFunctions.h>
#include <itkTimeProbe.h>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);
	
    TCLAP::ValueArg<std::string> t2Arg("l","t2","List of T2 relaxometry images",true,"","T2 relaxometry images",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask",false,"","computation mask",cmd);

    TCLAP::ValueArg<std::string> t1MapArg("","t1","T1 map",false,"","T1 map",cmd);
    TCLAP::ValueArg<std::string> b1MapArg("","b1","B1 map",false,"","B1 map",cmd);
    TCLAP::ValueArg<std::string> initialT2Arg("","init-t2","T2 initial map",false,"","T2 initial map",cmd);
    TCLAP::ValueArg<std::string> initialM0Arg("","init-m0","M0 initial map",false,"","M0 initial map",cmd);

    TCLAP::ValueArg<std::string> resT2Arg("o","out-t2","Result T2 image",true,"","result T2 image",cmd);
    TCLAP::ValueArg<std::string> resM0Arg("O","out-m0","Result M0 image",false,"","result M0 image",cmd);
	
    TCLAP::ValueArg<double> upperBoundT2Arg("u","upper-bound-t2","T2 value upper bound (default: 1000)",false,1000,"T2 value upper bound",cmd);
    TCLAP::ValueArg<double> upperBoundM0Arg("","upper-bound-m0","M0 value upper bound (default: 5000)",false,5000,"M0 value upper bound",cmd);

    TCLAP::ValueArg<double> echoSpacingArg("e","echo-spacing","Spacing between two successive echoes (default: 10)",false,10,"Spacing between echoes",cmd);
    TCLAP::ValueArg<double> excitationT2FlipAngleArg("","t2-ex-flip","Excitation flip angle for T2 (in degrees, default: 90)",false,90,"T2 excitation flip angle",cmd);
    TCLAP::ValueArg<double> t2FlipAngleArg("","t2-flip","All flip angles for T2 (in degrees, default: 180)",false,180,"T2 flip angle",cmd);
    TCLAP::ValueArg<double> backgroundSignalThresholdArg("t","signal-thr","Background signal threshold (default: 10)",false,10,"Background signal threshold",cmd);
    
    TCLAP::SwitchArg b1OnExcAngleArg("B", "b1-exc", "B1 is also applied to excitation angle",cmd,false);

    TCLAP::ValueArg<unsigned int> nbpArg("T","numberofthreads","Number of threads to run on (default : all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
	
    TCLAP::ValueArg<unsigned int> numOptimizerIterArg("","opt-iter","Maximal number of optimizer iterations (default: 200)",false,200,"Maximal number of optimizer iterations",cmd);
    TCLAP::ValueArg<double> optimizerStopConditionArg("","opt-stop","Optimizer stopping threshold (default: 1.0e-4)",false,1.0e-4,"Optimizer stopping threshold",cmd);
    TCLAP::ValueArg<double> optimizerInitialStepArg("i","opt-init","Optimizer initial step (default: 10)",false,10,"Optimizer initial step",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }
    
    typedef itk::Image <double,3> InputImageType;
    typedef itk::Image <double,4> Image4DType;
    typedef InputImageType OutputImageType;
    typedef anima::T2EPGRelaxometryEstimationImageFilter <InputImageType,OutputImageType> FilterType;
    typedef itk::ImageFileReader <InputImageType> InputImageReaderType;
    typedef itk::ImageFileWriter <OutputImageType> OutputImageWriterType;
    
    FilterType::Pointer mainFilter = FilterType::New();
	
    // Read and set T2 relaxometry images
    unsigned int numInputs = anima::setMultipleImageFilterInputsFromFileName<InputImageType,Image4DType,FilterType>(t2Arg.getValue(),mainFilter);
   
    mainFilter->SetEchoSpacing(echoSpacingArg.getValue());
    mainFilter->SetT2FlipAngles(t2FlipAngleArg.getValue() * M_PI / 180.0,numInputs);
    mainFilter->SetT2ExcitationFlipAngle(excitationT2FlipAngleArg.getValue() * M_PI / 180.0);

    mainFilter->SetM0UpperBound(upperBoundM0Arg.getValue());
    mainFilter->SetT2UpperBound(upperBoundT2Arg.getValue());
    
    mainFilter->SetMaximumOptimizerIterations(numOptimizerIterArg.getValue());
    mainFilter->SetOptimizerStopCondition(optimizerStopConditionArg.getValue());
    mainFilter->SetOptimizerInitialStep(optimizerInitialStepArg.getValue());

    if (t1MapArg.getValue() != "")
    {
        InputImageReaderType::Pointer t1MapRead = InputImageReaderType::New();
        t1MapRead->SetFileName(t1MapArg.getValue().c_str());
        t1MapRead->Update();
        
        mainFilter->SetT1Map(t1MapRead->GetOutput());
    }
    
    if (b1MapArg.getValue() != "")
    {
        InputImageReaderType::Pointer b1MapRead = InputImageReaderType::New();
        b1MapRead->SetFileName(b1MapArg.getValue().c_str());
        b1MapRead->Update();
        
        mainFilter->SetB1Map(b1MapRead->GetOutput());
    }
    
    if (initialT2Arg.getValue() != "")
    {
        InputImageReaderType::Pointer initT2Read = InputImageReaderType::New();
        initT2Read->SetFileName(initialT2Arg.getValue().c_str());
        initT2Read->Update();
        
        mainFilter->SetInitialT2Map(initT2Read->GetOutput());
    }
    
    if (initialM0Arg.getValue() != "")
    {
        InputImageReaderType::Pointer initM0Read = InputImageReaderType::New();
        initM0Read->SetFileName(initialM0Arg.getValue().c_str());
        initM0Read->Update();
        
        mainFilter->SetInitialM0Map(initM0Read->GetOutput());
    }

    if (maskArg.getValue() != "")
    {
        typedef itk::ImageFileReader < itk::Image <unsigned char, 3> > itkMaskReader;
        itkMaskReader::Pointer maskRead = itkMaskReader::New();
        maskRead->SetFileName(maskArg.getValue().c_str());
        maskRead->Update();
        
        mainFilter->SetComputationMask(maskRead->GetOutput());
    }
    
    mainFilter->SetB1OnExcitationAngle(b1OnExcAngleArg.isSet());

    mainFilter->SetAverageSignalThreshold(backgroundSignalThresholdArg.getValue());
    mainFilter->SetNumberOfThreads(nbpArg.getValue());
    
    itk::TimeProbe tmpTime;
    tmpTime.Start();
    
    mainFilter->Update();
    
    tmpTime.Stop();
    
    std::cout << "Total computation time: " << tmpTime.GetTotal() << std::endl;

    OutputImageWriterType::Pointer resultT2Writer = OutputImageWriterType::New();
    resultT2Writer->SetFileName(resT2Arg.getValue().c_str());
    resultT2Writer->SetUseCompression(true);
    resultT2Writer->SetInput(mainFilter->GetOutput(0));
    
    resultT2Writer->Update();
    
    if (resM0Arg.getValue() != "")
    {
        OutputImageWriterType::Pointer resultM0Writer = OutputImageWriterType::New();
        resultM0Writer->SetFileName(resM0Arg.getValue().c_str());
        resultM0Writer->SetUseCompression(true);
        resultM0Writer->SetInput(mainFilter->GetOutput(1));
        
        resultM0Writer->Update();
    }

    return 0;
}
