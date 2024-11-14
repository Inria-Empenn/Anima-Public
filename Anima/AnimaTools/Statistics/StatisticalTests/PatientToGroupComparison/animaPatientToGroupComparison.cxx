#include <iostream>
#include <tclap/CmdLine.h>

#include <animaReadWriteFunctions.h>
#include <animaPatientToGroupComparisonImageFilter.h>

//Update progression of the process
void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
    std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
}

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> refLTArg("i","input","Test Image",true,"","test image",cmd);
    TCLAP::ValueArg<std::string> dataLTArg("I","database","Database Image List",true,"","database image list",cmd);
    
    TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask",true,"","computation mask",cmd);
    TCLAP::ValueArg<std::string> resArg("o","outputname","Z-Score output image",true,"","Z-Score output image",cmd);
    TCLAP::ValueArg<std::string> resPValArg("O","outpvalname","P-value output image",true,"","P-Value output image",cmd);

    TCLAP::ValueArg<std::string> statTestArg("t","stat-test","Statistical test to use ([fisher],chi)",false,"fisher","statistical test",cmd);
    TCLAP::ValueArg<double> expVarArg("e","expvar","PCA threshold: threshold on eigenvalues to compute the new basis (default: 0.5)",false,0.5,"PCA threshold",cmd);
    TCLAP::ValueArg<unsigned int> numEigenArg("E","numeigenpca","Number of eigenvalues to keep (default: 6)",false,6,"Number of PCA eigen values",cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default: all cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }
    
    typedef itk::VectorImage<double,3> LogTensorImageType;
    typedef anima::PatientToGroupComparisonImageFilter<double> MAZScoreImageFilterType;

    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);

    MAZScoreImageFilterType::Pointer mainFilter = MAZScoreImageFilterType::New();
    mainFilter->SetComputationMask(anima::readImage < itk::Image <unsigned char, 3> > (maskArg.getValue()));
    mainFilter->SetNumberOfWorkUnits(nbpArg.getValue());

    mainFilter->SetInput(anima::readImage <LogTensorImageType> (refLTArg.getValue()));
    mainFilter->SetExplainedRatio(expVarArg.getValue());
    mainFilter->SetNumEigenValuesPCA(numEigenArg.getValue());

    if (statTestArg.getValue() == "chi")
        mainFilter->SetStatisticalTestType(MAZScoreImageFilterType::CHI_SQUARE);
    else
        mainFilter->SetStatisticalTestType(MAZScoreImageFilterType::FISHER);

    std::ifstream fileIn(dataLTArg.getValue());
    if (!fileIn.is_open())
    {
        std::cerr << "Could not open data file (" << dataLTArg.getValue() << ")" << std::endl;
        return EXIT_FAILURE;
    }
    
    while (!fileIn.eof())
    {
        char tmpStr[2048];
        fileIn.getline(tmpStr,2048);
        
        if (strcmp(tmpStr,"") == 0)
            continue;
        
        std::cout << "Loading tensor image " << tmpStr << "..." << std::endl;
        mainFilter->AddDatabaseInput(anima::readImage <LogTensorImageType> (tmpStr));
    }
    fileIn.close();

    mainFilter->AddObserver(itk::ProgressEvent(), callback);

    try
    {
        mainFilter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }
    
    anima::writeImage < itk::Image <double, 3> > (resArg.getValue(),mainFilter->GetOutput(0));
    anima::writeImage < itk::Image <double, 3> > (resPValArg.getValue(),mainFilter->GetOutput(1));

    return EXIT_SUCCESS;
}
