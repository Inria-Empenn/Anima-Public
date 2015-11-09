#include <iostream>
#include <tclap/CmdLine.h>

#include <animaLowMemNLMeansPatientToGroupComparisonBridge.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> refLTArg("i","input","Test Image",true,"","test image",cmd);
    TCLAP::ValueArg<std::string> dataLTArg("I","database","Database Image List",true,"","database image list",cmd);
    
	TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask",true,"","computation mask",cmd);
    TCLAP::ValueArg<std::string> resArg("o","outputpval","P-values output image",true,"","P-values output image",cmd);
    
    TCLAP::ValueArg<std::string> dbMeanDistAveArg("d","dbmeanave","DB average mean distance image",true,"","DB average mean distance image",cmd);
    TCLAP::ValueArg<std::string> dbMeanDistStdArg("D","dbmeanstd","DB std mean distance image",true,"","DB std mean distance image",cmd);
    TCLAP::ValueArg<std::string> dbCovDistAveArg("","dbcovave","DB average covariance distance image",true,"","DB average covariance distance image",cmd);
    TCLAP::ValueArg<std::string> dbCovDistStdArg("","dbcovstd","DB std covariance distance image",true,"","DB std covariance distance image",cmd);

    TCLAP::ValueArg<std::string> resScoreArg("O","outputscore","Score output image",false,"","Score output image",cmd);
    TCLAP::ValueArg<std::string> resNumPatchesArg("","outputnpatches","Number of patches output image",false,"","Number of patches output image",cmd);
	
	TCLAP::ValueArg<double> weightThrArg("w","weightthr","NL weight threshold: patches around have to be similar enough (default: 0.0)",false,0.0,"NL weight threshold",cmd);
	TCLAP::ValueArg<double> meanThrArg("M","patchmeanthr","Proportion of patches kept after mean test (default: 10.0)",false,10.0,"NL mean patch proportion",cmd);
	TCLAP::ValueArg<double> varThrArg("c","patchcovthr","Proportion of patches kept after covariance test (default: 10.0)",false,10.0,"NL covariance patch proportion",cmd);
	TCLAP::ValueArg<double> betaArg("b","beta","Beta parameter for local noise estimation (default: 1)",false,1,"Beta for local noise estimation",cmd);
    
    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default: all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
    
    TCLAP::ValueArg<unsigned int> patchHSArg("","patchhalfsize","Patch half size in each direction (default: 1)",false,1,"patch half size",cmd);
    TCLAP::ValueArg<unsigned int> patchSSArg("","patchstepsize","Patch step size for searching (default: 2)",false,2,"patch search step size",cmd);
    TCLAP::ValueArg<unsigned int> patchNeighArg("n","patchneighborhood","Patch half neighborhood size (default: 4)",false,4,"patch search neighborhood size",cmd);

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
	
	typedef anima::LowMemoryNLMeansPatientToGroupComparisonBridge NLMeansPatientToGroupComparisonBridgeType;
    
    NLMeansPatientToGroupComparisonBridgeType *mainFilter = new NLMeansPatientToGroupComparisonBridgeType;
    mainFilter->SetComputationMask(maskName);
    mainFilter->SetNumberOfThreads(nbpArg.getValue());
    
    mainFilter->SetPatchHalfSize(patchHSArg.getValue());
    mainFilter->SetSearchNeighborhood(patchNeighArg.getValue());
    mainFilter->SetSearchStepSize(patchSSArg.getValue());
    
    mainFilter->SetWeightThreshold(weightThrArg.getValue());
	mainFilter->SetMeanThreshold(meanThrArg.getValue());
	mainFilter->SetVarianceThreshold(varThrArg.getValue());
    
    mainFilter->SetBetaParameter(betaArg.getValue());
    
	mainFilter->SetDatabaseNames(dataName);
	mainFilter->SetTestFileName(refLTArg.getValue());

    mainFilter->SetDatabaseMeanDistanceAverageFileName(dbMeanDistAveArg.getValue());
	mainFilter->SetDatabaseMeanDistanceStdFileName(dbMeanDistStdArg.getValue());
	mainFilter->SetDatabaseCovarianceDistanceAverageFileName(dbCovDistAveArg.getValue());
	mainFilter->SetDatabaseCovarianceDistanceStdFileName(dbCovDistStdArg.getValue());
    
	mainFilter->SetOutputScoreName(resScoreArg.getValue());
	mainFilter->SetOutputPValName(resArg.getValue());
	mainFilter->SetOutputNPatchesName(resNumPatchesArg.getValue());
	
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
