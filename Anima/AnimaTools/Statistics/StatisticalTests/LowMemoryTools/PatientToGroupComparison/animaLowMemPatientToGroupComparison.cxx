#include <iostream>
#include <tclap/CmdLine.h>

#include "animaLowMemPatientToGroupComparisonBridge.h"

using namespace std;

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> refLTArg("i","inputlogtens","Log-Tensor Test Image",true,"","log-tensor test image",cmd);    
    TCLAP::ValueArg<std::string> dataLTArg("I","databaselogtens","Log-Tensor Database Image List",true,"","log-tensor database image list",cmd);
    
	TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask",true,"","computation mask",cmd);
    TCLAP::ValueArg<std::string> resArg("o","outputname","output image prefix",true,"","output image prefix",cmd);
    TCLAP::ValueArg<std::string> resPValArg("O","outpvalname","P-value output image",true,"","P-Value output image",cmd);

    TCLAP::ValueArg<std::string> statTestArg("t","stat-test","Statistical test to use ([fisher],chi)",false,"fisher","statistical test",cmd);
    TCLAP::ValueArg<double> expVarArg("e","expvar","PCA threshold: threshold on eigenvalues to compute the new basis (default: 0.9)",false,0.9,"PCA threshold",cmd);
    TCLAP::ValueArg<unsigned int> numEigenArg("E","numeigenpca","Number of eigenvalues to keep (default: 6)",false,6,"Number of PCA eigen values",cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default : all cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
    
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
    
    string dataLTName, maskName;
    dataLTName = dataLTArg.getValue();
    maskName = maskArg.getValue();
	
    typedef anima::LowMemoryPatientToGroupComparisonBridge MultiAtlasZScoreBridgeType;
    
    MultiAtlasZScoreBridgeType *mainFilter = new MultiAtlasZScoreBridgeType;
    mainFilter->SetComputationMask(maskName);
    mainFilter->SetNumberOfWorkUnits(nbpArg.getValue());

    mainFilter->SetNumEigenValuesPCA(numEigenArg.getValue());
    mainFilter->SetExplainedRatio(expVarArg.getValue());
	
	mainFilter->SetDataLTFileNames(dataLTName);
	mainFilter->SetTestLTFileName(refLTArg.getValue());
    
    mainFilter->SetOutputName(resArg.getValue());
    mainFilter->SetOutputPValName(resPValArg.getValue());

	mainFilter->SetNbSplits(splitsArg.getValue());

    if (statTestArg.getValue() == "chi")
        mainFilter->SetStatisticalTestType(MultiAtlasZScoreBridgeType::MainFilterType::CHI_SQUARE);
    else
        mainFilter->SetStatisticalTestType(MultiAtlasZScoreBridgeType::MainFilterType::FISHER);

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
