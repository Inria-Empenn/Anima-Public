#include <iostream>
#include <tclap/CmdLine.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <animaLowMemCramersTestBridge.h>
#include <itkTimeProbe.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
	
    TCLAP::ValueArg<std::string> logsArg("l","logfilelist","File containing the list of log-tensors",true,"","log-tensors list",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask",true,"","computation mask",cmd);
    TCLAP::ValueArg<std::string> outlierMasksArg("M","outliermasks","Outlier masks list",false,"","outlier masks list",cmd);
    TCLAP::ValueArg<std::string> resArg("o","outputprefix","Result image prefix",true,"","result image prefix",cmd);
	
    TCLAP::ValueArg<std::string> fGrpArg("","fg","Text file containing the indexes of images from the first group",true,"","first group indexes list",cmd);
    TCLAP::ValueArg<std::string> sGrpArg("","sg","Text file containing the indexes of images from the second group",true,"","second group indexes list",cmd);
	
    TCLAP::ValueArg<unsigned int> nbSamplesArg("n","nbbootstrapsamples","Number of permutation samples (default: 5000)",false,5000,"# permutation samples",cmd);
    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default : all cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
    
    TCLAP::ValueArg<unsigned int> splitsArg("s","split","Split image for low memory (default: 2)",false,2,"Number of splits",cmd);
    TCLAP::ValueArg<int> specSplitArg("S","splittoprocess","Specific split to process (use to run on cluster (default: -1 = all)",false,-1,"Split to process",cmd);
    TCLAP::SwitchArg genOutputDescroArg("G","generateouputdescription","Generate ouptut description data",cmd,false);;
	
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }
	
    typedef anima::LowMemoryCramersTestBridge CramersBridgeType;
    
    std::string logsName, maskName, resName;
    logsName = logsArg.getValue();
    maskName = maskArg.getValue();
    resName = resArg.getValue();
	
    unsigned int nbSamples = nbSamplesArg.getValue();
	
    itk::TimeProbe tmpTime;
    tmpTime.Start();
    
    std::string fGrpName, sGrpName;
    fGrpName = fGrpArg.getValue();
    sGrpName = sGrpArg.getValue();
    
    CramersBridgeType *cramersFilter = new CramersBridgeType;
    cramersFilter->SetComputationMask(maskName);
    cramersFilter->SetNbSamples(nbSamples);
        
    cramersFilter->SetInputFileNames(logsName);
    if (outlierMasksArg.getValue() != "")
        cramersFilter->SetOutlierMaskFileNames(outlierMasksArg.getValue());
    
    cramersFilter->SetNbSplits(splitsArg.getValue());
    cramersFilter->SetOutputPrefix(resName);
    cramersFilter->SetIndexesFromFiles(fGrpName,sGrpName);
	
    unsigned int nbProcs = nbpArg.getValue();
    cramersFilter->SetNumberOfThreads(nbProcs);
        
    cramersFilter->Update(specSplitArg.getValue(),genOutputDescroArg.getValue());
    
    tmpTime.Stop();
    
    std::cout << "Total computation time: " << tmpTime.GetTotal() << std::endl;
	
	delete cramersFilter;
	
    return 0;
}
