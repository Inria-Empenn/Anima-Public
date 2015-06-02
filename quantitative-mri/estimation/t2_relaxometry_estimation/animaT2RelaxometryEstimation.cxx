#include <iostream>
#include <tclap/CmdLine.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <animaT2RelaxometryEstimationImageFilter.h>
#include <animaReadWriteFunctions.h>
#include <itkTimeProbe.h>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',"1.0");
	
    TCLAP::ValueArg<std::string> t2Arg("l","t2","List of T2 relaxometry images",true,"","T2 relaxometry images",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask",false,"","computation mask",cmd);

    TCLAP::ValueArg<std::string> t1MapArg("","t1","T1 map",false,"","T1 map",cmd);

    TCLAP::ValueArg<std::string> resT2Arg("o","out-t2","Result T2 image",true,"","result T2 image",cmd);
    TCLAP::ValueArg<std::string> resM0Arg("O","out-m0","Result M0 image",false,"","result M0 image",cmd);
	
    TCLAP::ValueArg<double> trArg("","tr","Repetition time for T2 relaxometry (default: 5000)",false,5000,"repetition time",cmd);
    TCLAP::ValueArg<double> upperBoundT2Arg("u","upper-bound-t2","T2 value upper bound (default: 1000)",false,1000,"T2 value upper bound",cmd);
    TCLAP::ValueArg<double> upperBoundM0Arg("","upper-bound-m0","M0 value upper bound (default: 5000)",false,5000,"M0 value upper bound",cmd);

    TCLAP::ValueArg<std::string> echoSpacingArg("e","echo-spacing","Spacing between two successive echoes (default: 10), or a file containing the echo times for each acquisition",false,"10","Spacing between echoes",cmd);
    TCLAP::ValueArg<double> backgroundSignalThresholdArg("t","signal-thr","Background signal threshold (default: 10)",false,10,"Background signal threshold",cmd);
    
    TCLAP::ValueArg<unsigned int> nbpArg("T","numberofthreads","Number of threads to run on (default : all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
	
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
    typedef anima::T2RelaxometryEstimationImageFilter <InputImageType,OutputImageType> FilterType;
    typedef itk::ImageFileReader <InputImageType> InputImageReaderType;
    typedef itk::ImageFileWriter <OutputImageType> OutputImageWriterType;
    
    FilterType::Pointer mainFilter = FilterType::New();
	
    // Read and set T2 relaxometry images
    anima::setMultipleImageFilterInputsFromFileName<InputImageType,Image4DType,FilterType>(t2Arg.getValue(),mainFilter);
   
    mainFilter->SetTRValue(trArg.getValue());


    // get the echo time from double or from file
    unsigned int numInputs=mainFilter->GetNumberOfIndexedInputs();
    std::vector<double> echoTime(numInputs);
    try
    {
    const double echoSpacing=boost::lexical_cast<double>(echoSpacingArg.getValue());
    for(unsigned int index=0;index<echoTime.size();++index)
        echoTime[index]=(index+1)*(echoSpacing);
    }
    catch(const boost::bad_lexical_cast &)  // then try to read file
    {
        std::ifstream inputFile(echoSpacingArg.getValue().c_str());
        if (!inputFile.is_open())
        {
            std::cerr << "Entry for echo spacing is neither a double, neither a readable file" << echoSpacingArg.getValue() << std::endl;
            return EXIT_FAILURE;
        }

        unsigned int index = 0;
        while (!inputFile.eof())
        {
           inputFile>>echoTime[index];
           ++index;
        }
        if(index!=numInputs)
        {
            std::cerr << "Number of echos in " << echoSpacingArg.getValue() << "different from number of input acquisitions."<< std::endl;
            return EXIT_FAILURE;
        }
    }
    mainFilter->SetEchoTime(echoTime);
    //

    mainFilter->SetM0UpperBoundValue(upperBoundM0Arg.getValue());
    mainFilter->SetT2UpperBoundValue(upperBoundT2Arg.getValue());
    
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
