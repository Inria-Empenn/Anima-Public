#include <iostream>
#include <tclap/CmdLine.h>

#include <animaT1SERelaxometryEstimationImageFilter.h>
#include <animaReadWriteFunctions.h>
#include <itkTimeProbe.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);
	
    TCLAP::ValueArg<std::string> t1Arg("l","t1","List of T1 relaxometry images",true,"","T1 relaxometry images",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask",false,"","computation mask",cmd);
    
    TCLAP::ValueArg<std::string> resT1Arg("o","out-t1","Result T1 image",true,"","result T1 image",cmd);
    TCLAP::ValueArg<std::string> resM0Arg("O","out-m0","Result M0 image",false,"","result M0 image",cmd);
	
    TCLAP::ValueArg<double> upperBoundT1Arg("u","upper-bound-t1","T1 value upper bound (default: 5000)",false,5000,"T1 value upper bound",cmd);
    TCLAP::ValueArg<double> upperBoundM0Arg("","upper-bound-m0","M0 value upper bound (default: 5000)",false,5000,"M0 value upper bound",cmd);

    TCLAP::ValueArg<std::string> trArg("","tr","Text file listing repetition times for input images",true,"","repetition times list",cmd);
    TCLAP::ValueArg<double> backgroundSignalThresholdArg("t","signal-thr","Background signal threshold (default: 10)",false,10,"Background signal threshold",cmd);
    
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
        return EXIT_FAILURE;
    }
    
    typedef itk::Image <double,3> InputImageType;
    typedef InputImageType OutputImageType;
    typedef anima::T1SERelaxometryEstimationImageFilter <InputImageType,OutputImageType> FilterType;
    typedef itk::ImageFileWriter <OutputImageType> OutputImageWriterType;
    
    FilterType::Pointer mainFilter = FilterType::New();
	
    mainFilter->SetM0UpperBound(upperBoundM0Arg.getValue());
    mainFilter->SetT1UpperBound(upperBoundT1Arg.getValue());
    
    mainFilter->SetMaximumOptimizerIterations(numOptimizerIterArg.getValue());
    mainFilter->SetOptimizerStopCondition(optimizerStopConditionArg.getValue());
    mainFilter->SetOptimizerInitialStep(optimizerInitialStepArg.getValue());

    // Read and set T1 relaxometry images
    anima::setMultipleImageFilterInputsFromFileName<InputImageType,FilterType>(t1Arg.getValue(),mainFilter);
    
    std::ifstream inputTRValues(trArg.getValue().c_str());
    if (!inputTRValues.is_open())
    {
        std::cerr << "Text file for TR values not readable " << trArg.getValue() << std::endl;
        return EXIT_FAILURE;
    }

    std::vector <double> trValues;
    while (!inputTRValues.eof())
    {
        char tmpStr[8192];
        inputTRValues.getline(tmpStr,8192);
        std::string workStr(tmpStr);
        if (workStr == "")
            continue;
        
        workStr.erase(workStr.find_last_not_of(" \n\r\t")+1);
        std::istringstream iss(workStr);
        std::string shortStr;
        iss >> shortStr;
        trValues.push_back(std::stod(shortStr));
    }
    
    mainFilter->SetTRValues(trValues);
    
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
    
    itk::TimeProbe tmpTime;
    tmpTime.Start();
    
    mainFilter->Update();
    
    tmpTime.Stop();
    
    std::cout << "Total computation time: " << tmpTime.GetTotal() << std::endl;

    anima::writeImage <OutputImageType> (resT1Arg.getValue(),mainFilter->GetOutput(1));
    
    if (resM0Arg.getValue() != "")
        anima::writeImage <OutputImageType> (resM0Arg.getValue(),mainFilter->GetOutput(0));

    return EXIT_SUCCESS;
}
