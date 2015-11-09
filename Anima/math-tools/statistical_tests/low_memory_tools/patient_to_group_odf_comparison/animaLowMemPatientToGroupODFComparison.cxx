#include <iostream>
#include <tclap/CmdLine.h>

#include <animaLowMemPatientToGroupODFComparisonBridge.h>
#include <animaODFFunctions.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);
	
    TCLAP::ValueArg<std::string> refODFArg("i","inputodf","ODF Test Image",true,"","ODF test image",cmd);
    TCLAP::ValueArg<std::string> dataODFArg("I","databaseodf","ODF Database Image List",true,"","ODF database image list",cmd);
    
	TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask",true,"","computation mask",cmd);
    TCLAP::ValueArg<std::string> resArg("o","outputname","Z-Score output image",true,"","Z-Score output image",cmd);
    TCLAP::ValueArg<std::string> resPValArg("O","outpvalname","P-value output image",true,"","P-Value output image",cmd);

    TCLAP::ValueArg<std::string> statTestArg("t","stat-test","Statistical test to use ([fisher],chi)",false,"fisher","statistical test",cmd);
    TCLAP::ValueArg<double> expVarArg("e","expvar","PCA threshold: threshold on eigenvalues to compute the new basis (default: 0.9)",false,0.9,"PCA threshold",cmd);
	TCLAP::ValueArg<unsigned int> numEigenArg("E","numeigenpca","Number of eigenvalues to keep (default: 6)",false,6,"Number of PCA eigen values",cmd);
    
    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default: all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);	
	TCLAP::ValueArg<unsigned int> nbThetaArg("T","theta","Number of theta values (theta varies between 0 and pi/2",false,0,"number of theta values",cmd);
    TCLAP::ValueArg<unsigned int> nbPhiArg("P","phi","Number of phi values (theta varies between 0 and 2 pi",false,0,"number of phi values",cmd);
	
	TCLAP::ValueArg<std::string> samplesFileNameArg("d","sampledirectionsfile","Samples directions in a text file",false,"","Samples directions file",cmd);
	
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
	
    std::string dataODFName, maskName;
    dataODFName = dataODFArg.getValue();
    maskName = maskArg.getValue();
	
    typedef anima::LowMemoryPatientToGroupODFComparisonBridge MultiAtlasZScoreBridgeType;
    
    MultiAtlasZScoreBridgeType *mainFilter = new MultiAtlasZScoreBridgeType;
    mainFilter->SetComputationMask(maskName);
    mainFilter->SetNumberOfThreads(nbpArg.getValue());
	
	mainFilter->SetNumEigenValuesPCA(numEigenArg.getValue());
	mainFilter->SetExplainedRatio(expVarArg.getValue());
	
	mainFilter->SetDataODFFileNames(dataODFName);
	mainFilter->SetTestODFFileName(refODFArg.getValue());
    
    std::vector < std::vector <double> > sampleDirections = anima::InitializeSampleDirections(nbThetaArg.getValue(),nbPhiArg.getValue(),
                                                                                              samplesFileNameArg.getValue());
	
	mainFilter->SetSampleDirections(sampleDirections);

    if (statTestArg.getValue() == "chi")
        mainFilter->SetStatisticalTestType(MultiAtlasZScoreBridgeType::MainFilterType::CHI_SQUARE);
    else
        mainFilter->SetStatisticalTestType(MultiAtlasZScoreBridgeType::MainFilterType::FISHER);

	mainFilter->SetOutputName(resArg.getValue());
    mainFilter->SetOutputPValName(resPValArg.getValue());

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

