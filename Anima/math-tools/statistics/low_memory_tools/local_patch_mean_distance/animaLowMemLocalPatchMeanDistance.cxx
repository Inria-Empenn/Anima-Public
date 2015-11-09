#include <iostream>
#include <tclap/CmdLine.h>

#include <animaLowMemLocalPatchMeanDistanceBridge.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> dataLTArg("i","database","Database Image List",true,"","database image list",cmd);
    
	TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask",true,"","computation mask",cmd);
    TCLAP::ValueArg<std::string> resArg("o","outputmean","Average distance output image",true,"","Average distance output image",cmd);
    
    TCLAP::ValueArg<std::string> resStdArg("O","outputstd","Standard deviation output image",false,"","Standard deviation output image",cmd);
    
    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default: all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
    
    TCLAP::ValueArg<unsigned int> patchHSArg("","patchhalfsize","Patch half size in each direction (default: 1)",false,1,"patch half size",cmd);

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
        return(1);
    }
    
    std::string dataName, maskName;
    dataName = dataLTArg.getValue();
    maskName = maskArg.getValue();
	
	typedef anima::LowMemoryLocalPatchMeanDistanceBridge MainBridgeType;
    
    MainBridgeType *mainFilter = new MainBridgeType;
    mainFilter->SetComputationMask(maskName);
    mainFilter->SetNumberOfThreads(nbpArg.getValue());
    
    mainFilter->SetPatchHalfSize(patchHSArg.getValue());
	mainFilter->SetDatabaseNames(dataName);
    
	mainFilter->SetOutputMeanName(resArg.getValue());
	mainFilter->SetOutputStdName(resStdArg.getValue());
	
	mainFilter->SetNbSplits(splitsArg.getValue());

    try
    {
        mainFilter->Update(specSplitArg.getValue(),genOutputDescroArg.getValue());
    }
    catch(itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return 1;
    }
    
	delete mainFilter;
	
    return 0;
}
