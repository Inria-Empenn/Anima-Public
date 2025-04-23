#include <iostream>
#include <tclap/CmdLine.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <animaCramersTestImageFilter.h>
#include <itkTimeProbe.h>
#include <itkCommand.h>

//Update progression of the process
void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
    std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
}

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
	
    TCLAP::ValueArg<std::string> logsArg("l","logfilelist","File containing the list of log-tensors",true,"","log-tensors list",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask",true,"","computation mask",cmd);
    TCLAP::ValueArg<std::string> outlierMasksArg("M","outliermasks","Outlier masks list",false,"","outlier masks list",cmd);
    TCLAP::ValueArg<std::string> resArg("o","outputname","Result image",true,"","result image",cmd);
	
    TCLAP::ValueArg<std::string> fGrpArg("f","firstgroupfile","Text file containing the indexes of images from the first group",true,"","first group indexes list",cmd);
    TCLAP::ValueArg<std::string> sGrpArg("s","secondgroupfile","Text file containing the indexes of images from the second group",true,"","second group indexes list",cmd);
	
    TCLAP::ValueArg<unsigned int> nbSamplesArg("n","nbbootstrapsamples","Number of permutation samples (default: 5000)",false,5000,"# permutation samples",cmd);
    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default : all cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
	
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }
	
    std::string logsName, maskName, resName;
    logsName = logsArg.getValue();
    maskName = maskArg.getValue();
    resName = resArg.getValue();
	
    unsigned int nbSamples = nbSamplesArg.getValue();
	
    itk::TimeProbe tmpTime;
    tmpTime.Start();
	
    typedef itk::ImageFileWriter < itk::Image <double, 3> > itkOutputWriter;
    typedef itk::ImageFileReader < itk::Image <unsigned char, 3> > itkMaskReader;
    itkMaskReader::Pointer maskRead = itkMaskReader::New();
    maskRead->SetFileName(maskName.c_str());
    maskRead->Update();
	
    std::vector <unsigned int> firstGroupIndexes, secGroupIndexes;
    std::string fGrpName, sGrpName;
    fGrpName = fGrpArg.getValue();
    sGrpName = sGrpArg.getValue();
	
    if (strcmp(fGrpName.c_str(),"") != 0)
    {
        std::ifstream fGrpfile(fGrpName.c_str());
        while (!fGrpfile.eof())
        {
            unsigned int tmpVal;
            char tmpStr[2048];
            fGrpfile.getline(tmpStr,2048);
            
            if (strcmp(tmpStr,"") != 0)
            {
				sscanf(tmpStr,"%d",&tmpVal);
				firstGroupIndexes.push_back(tmpVal);
            }
        }
		
        fGrpfile.close();
    }
	
    if (strcmp(sGrpName.c_str(),"") != 0)
    {
        std::ifstream sGrpfile(sGrpName.c_str());
        while (!sGrpfile.eof())
        {
            unsigned int tmpVal;
            char tmpStr[2048];
            sGrpfile.getline(tmpStr,2048);
            
            if (strcmp(tmpStr,"") != 0)
            {
				sscanf(tmpStr,"%d",&tmpVal);
				secGroupIndexes.push_back(tmpVal);
            }
        }
		
        sGrpfile.close();
    }
	
    typedef itk::VectorImage<double,3> LogTensorImageType;
    typedef itk::ImageFileReader <LogTensorImageType> itkInputReader;
    typedef anima::CramersTestImageFilter<double> CramersFilterType;
    
    CramersFilterType::Pointer cramersFilter = CramersFilterType::New();
    cramersFilter->SetComputationMask(maskRead->GetOutput());
    cramersFilter->SetNbSamples(nbSamples);
        
    std::ifstream fileIn(logsName.c_str());
    
    int nbPats = 0;
    itkInputReader::Pointer imageReader;
    
    while (!fileIn.eof())
    {
        char tmpStr[2048];
        fileIn.getline(tmpStr,2048);
        
        if (strcmp(tmpStr,"") == 0)
            continue;
        
        std::cout << "Loading image " << nbPats << " " << tmpStr << "..." << std::endl;
        imageReader = itkInputReader::New();
        imageReader->SetFileName(tmpStr);
        imageReader->Update();
        
        cramersFilter->SetInput(nbPats,imageReader->GetOutput());        
		nbPats++;
    }
    fileIn.close();
    
    if (outlierMasksArg.getValue() != "")
    {
        std::ifstream masksFile(outlierMasksArg.getValue().c_str());
        itkMaskReader::Pointer mReader;
        
        while (!masksFile.eof())
        {
            char tmpStr[2048];
            masksFile.getline(tmpStr,2048);
            
            if (strcmp(tmpStr,"") == 0)
                continue;
            
            mReader = itkMaskReader::New();
            mReader->SetFileName(tmpStr);
            mReader->Update();
            
            cramersFilter->AddOutlierMask(mReader->GetOutput());
        }
        
        masksFile.close();
    }
    
    cramersFilter->SetFirstGroupIndexes(firstGroupIndexes);
    cramersFilter->SetSecondGroupIndexes(secGroupIndexes);
	
    unsigned int nbProcs = nbpArg.getValue();
    cramersFilter->SetNumberOfWorkUnits(nbProcs);

    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);
    cramersFilter->AddObserver(itk::ProgressEvent(), callback );

    cramersFilter->Update();
    
    itkOutputWriter::Pointer resultWriter = itkOutputWriter::New();
    resultWriter->SetFileName(resName.c_str());
    resultWriter->SetUseCompression(true);
    resultWriter->SetInput(cramersFilter->GetOutput());
    
    resultWriter->Update();
    
    tmpTime.Stop();
    
    std::cout << "Total computation time: " << tmpTime.GetTotal() << std::endl;
	
    return 0;
}
