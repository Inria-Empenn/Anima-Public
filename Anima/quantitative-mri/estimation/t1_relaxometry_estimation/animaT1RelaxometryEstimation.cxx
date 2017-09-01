#include <iostream>
#include <tclap/CmdLine.h>

#include <cmath>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <animaT1RelaxometryEstimationImageFilter.h>
#include <animaReadWriteFunctions.h>
#include <itkTimeProbe.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);
	
    TCLAP::ValueArg<std::string> t1Arg("l","t1","List of T1 relaxometry images",true,"","T1 relaxometry images",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask",false,"","computation mask",cmd);
    TCLAP::ValueArg<std::string> b1Arg("b","b1","B1 inhomogenity map",false,"","B1 inhomogenity map",cmd);
    
    TCLAP::ValueArg<std::string> resT1Arg("o","out-t1","Result T1 image",true,"","result T1 image",cmd);
    TCLAP::ValueArg<std::string> resM0Arg("O","out-m0","Result M0 image",false,"","result M0 image",cmd);
	
    TCLAP::ValueArg<double> trArg("","tr","Repetition time for T1 relaxometry (default: 15)",false,15,"repetition time",cmd);
    TCLAP::ValueArg<std::string> flipAnglesArg("f","flip","Text file listing flip angles (in degrees)",true,"","T1 flip angles list",cmd);
    
    TCLAP::ValueArg<double> upperBoundT1Arg("u","upper-bound-t1","T1 value upper bound (default: 5000)",false,5000,"T1 upper bound",cmd);
    TCLAP::ValueArg<double> upperBoundM0Arg("","upper-bound-m0","M0 value upper bound (default: 5000)",false,5000,"M0 upper bound",cmd);
    TCLAP::ValueArg<double> backgroundSignalThresholdArg("t","signal-thr","Background signal threshold (default: 10)",false,10,"Background signal threshold",cmd);
    
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
    typedef anima::T1RelaxometryEstimationImageFilter <InputImageType,OutputImageType> FilterType;
    typedef itk::ImageFileReader <InputImageType> InputImageReaderType;
    typedef itk::ImageFileWriter <OutputImageType> OutputImageWriterType;
    
    FilterType::Pointer mainFilter = FilterType::New();
	
    // Read and set T1 relaxometry images    
    anima::setMultipleImageFilterInputsFromFileName<InputImageType,FilterType>(t1Arg.getValue(),mainFilter);
    
    std::ifstream inputFlip(flipAnglesArg.getValue().c_str());
    if (!inputFlip.is_open())
    {
        std::cerr << "Text file for flip angles not readable " << flipAnglesArg.getValue() << std::endl;
        return EXIT_FAILURE;
    }

    std::vector <double> flipAngles;
    while (!inputFlip.eof())
    {
        char tmpStr[8192];
        inputFlip.getline(tmpStr,8192);
        std::string workStr(tmpStr);
        if (workStr == "")
            continue;
        
        workStr.erase(workStr.find_last_not_of(" \n\r\t")+1);
        std::istringstream iss(workStr);
        std::string shortStr;
        iss >> shortStr;
        flipAngles.push_back(std::stod(shortStr) * M_PI / 180.0);
    }
    
    mainFilter->SetFlipAngles(flipAngles);
    mainFilter->SetTRValue(trArg.getValue());
    mainFilter->SetM0UpperBoundValue(upperBoundM0Arg.getValue());
    mainFilter->SetT1UpperBoundValue(upperBoundT1Arg.getValue());
    
    if (b1Arg.getValue() != "")
    {
        InputImageReaderType::Pointer b1Read = InputImageReaderType::New();
        b1Read->SetFileName(b1Arg.getValue().c_str());
        b1Read->Update();
        
        mainFilter->SetB1Map(b1Read->GetOutput());
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
    
    itk::TimeProbe tmpTime;
    tmpTime.Start();
    
    mainFilter->Update();
    
    tmpTime.Stop();
    
    std::cout << "Total computation time: " << tmpTime.GetTotal() << std::endl;

    OutputImageWriterType::Pointer resultT1Writer = OutputImageWriterType::New();
    resultT1Writer->SetFileName(resT1Arg.getValue().c_str());
    resultT1Writer->SetUseCompression(true);
    resultT1Writer->SetInput(mainFilter->GetOutput(0));
    
    resultT1Writer->Update();
    
    if (resM0Arg.getValue() != "")
    {
        OutputImageWriterType::Pointer resultM0Writer = OutputImageWriterType::New();
        resultM0Writer->SetFileName(resM0Arg.getValue().c_str());
        resultM0Writer->SetUseCompression(true);
        resultM0Writer->SetInput(mainFilter->GetOutput(1));
        
        resultM0Writer->Update();
    }

    return EXIT_SUCCESS;
}
