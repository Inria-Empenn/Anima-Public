#include <iostream>
#include <tclap/CmdLine.h>

#include <animaReadWriteFunctions.h>
#include <animaLocalPatchCovarianceDistanceImageFilter.h>
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
    
    TCLAP::ValueArg<std::string> dataLTArg("i","database","Database Image List",true,"","database image list",cmd);
    
	TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask",true,"","computation mask",cmd);
    TCLAP::ValueArg<std::string> resArg("o","outputmean","Average distance output image",true,"","Average distance output image",cmd);
    
    TCLAP::ValueArg<std::string> resStdArg("O","outputstd","Standard deviation output image",false,"","Standard deviation output image",cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default: all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    TCLAP::ValueArg<unsigned int> patchHSArg("","patchhalfsize","Patch half size in each direction (default: 1)",false,1,"patch half size",cmd);
    
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
    typedef anima::LocalPatchCovarianceDistanceImageFilter<double> MainFilterType;
       
    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);

    MainFilterType::Pointer mainFilter = MainFilterType::New();
    mainFilter->SetComputationMask(anima::readImage < itk::Image <unsigned char, 3> > (maskArg.getValue()));
    mainFilter->SetNumberOfThreads(nbpArg.getValue());

	mainFilter->SetPatchHalfSize(patchHSArg.getValue());
    mainFilter->AddObserver(itk::ProgressEvent(), callback);

    std::ifstream fileIn(dataLTArg.getValue());
    
    unsigned int count = 0;
    while (!fileIn.eof())
    {
        char tmpStr[2048];
        fileIn.getline(tmpStr,2048);
        
        if (strcmp(tmpStr,"") == 0)
            continue;
        
        std::cout << "Loading database image " << tmpStr << "..." << std::endl;
        mainFilter->SetInput(count,anima::readImage <LogTensorImageType> (tmpStr));
        ++count;
    }
    fileIn.close();
  	
    try
    {
        mainFilter->Update();
    }
    catch(itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }
    
    anima::writeImage < itk::Image <double, 3> > (resArg.getValue(),mainFilter->GetOutput(0));
    
    if (resStdArg.getValue() != "")
        anima::writeImage < itk::Image <double, 3> > (resStdArg.getValue(),mainFilter->GetOutput(1));
    
    return EXIT_SUCCESS;
}
