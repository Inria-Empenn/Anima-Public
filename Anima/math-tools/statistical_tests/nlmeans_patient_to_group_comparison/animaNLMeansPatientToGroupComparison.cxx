#include <iostream>
#include <tclap/CmdLine.h>

#include <animaReadWriteFunctions.h>
#include <animaNLMeansPatientToGroupComparisonImageFilter.h>
#include <itkTimeProbe.h>

//Update progression of the process
void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
    std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
}

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
    TCLAP::ValueArg<double> meanThrArg("M","patchmeanthr","Tolerance for means test (test if meansTest > meanDatabase + M * stdDatabase, default: M=2.5)",false,2.5,"NL mean patch proportion",cmd);
    TCLAP::ValueArg<double> varThrArg("c","patchcovthr","Tolerance for covariance test (test if covDist > meanDatabase + c * stdDatabase, default: c=2.5)",false,2.5,"NL covariance patch proportion",cmd);

    TCLAP::ValueArg<double> betaArg("b","beta","Beta parameter for local noise estimation (default: 1)",false,1,"Beta for local noise estimation",cmd);
    
    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default: all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    TCLAP::ValueArg<unsigned int> patchHSArg("","patchhalfsize","Patch half size in each direction (default: 1)",false,1,"patch half size",cmd);
    TCLAP::ValueArg<unsigned int> patchSSArg("","patchstepsize","Patch step size for searching (default: 1)",false,1,"patch search step size",cmd);
    TCLAP::ValueArg<unsigned int> patchNeighArg("n","patchneighborhood","Patch half neighborhood size (default: 2)",false,2,"patch search neighborhood size",cmd);
    
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
    typedef anima::NLMeansPatientToGroupComparisonImageFilter<double> NLComparisonImageFilterType;
    
    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);

    NLComparisonImageFilterType::Pointer mainFilter = NLComparisonImageFilterType::New();
    mainFilter->SetComputationMask(anima::readImage < itk::Image <unsigned char, 3> > (maskArg.getValue()));
    mainFilter->SetNumberOfThreads(nbpArg.getValue());

    mainFilter->SetWeightThreshold(weightThrArg.getValue());
    mainFilter->SetMeanThreshold(meanThrArg.getValue());
    mainFilter->SetVarianceThreshold(varThrArg.getValue());

    mainFilter->SetPatchHalfSize(patchHSArg.getValue());
    mainFilter->SetSearchStepSize(patchSSArg.getValue());
    mainFilter->SetSearchNeighborhood(patchNeighArg.getValue());
    mainFilter->SetBetaParameter(betaArg.getValue());

    mainFilter->SetInput(0,anima::readImage <LogTensorImageType> (refLTArg.getValue()));
    mainFilter->AddObserver(itk::ProgressEvent(), callback);

    std::ifstream fileIn(dataLTArg.getValue());
    
    while (!fileIn.eof())
    {
        char tmpStr[2048];
        fileIn.getline(tmpStr,2048);
        
        if (strcmp(tmpStr,"") == 0)
            continue;
        
        std::cout << "Loading database image " << tmpStr << "..." << std::endl;
        mainFilter->AddDatabaseInput(anima::readImage <LogTensorImageType> (tmpStr));
    }
    fileIn.close();
    
    mainFilter->SetDatabaseMeanDistanceAverage(anima::readImage < itk::Image <double, 3> > (dbMeanDistAveArg.getValue()));
    mainFilter->SetDatabaseMeanDistanceStd(anima::readImage < itk::Image <double, 3> > (dbMeanDistStdArg.getValue()));

    mainFilter->SetDatabaseCovarianceDistanceAverage(anima::readImage < itk::Image <double, 3> > (dbCovDistAveArg.getValue()));
    mainFilter->SetDatabaseCovarianceDistanceStd(anima::readImage < itk::Image <double, 3> > (dbCovDistStdArg.getValue()));

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
    
    if (resScoreArg.getValue() != "")
        anima::writeImage < itk::Image <double, 3> > (resScoreArg.getValue(),mainFilter->GetOutput(1));

    if (resNumPatchesArg.getValue() != "")
        anima::writeImage < itk::Image <double, 3> > (resNumPatchesArg.getValue(),mainFilter->GetOutput(2));
    
    return EXIT_SUCCESS;
}
