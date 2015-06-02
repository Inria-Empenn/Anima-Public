#include <iostream>
#include <tclap/CmdLine.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <animaLocalPatchCovarianceDistanceImageFilter.h>
#include <itkTimeProbe.h>

using namespace std;

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',"1.0");
    
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
        return(1);
    }
    
    string dataLTName, maskName;
    dataLTName = dataLTArg.getValue();
    maskName = maskArg.getValue();
    
    typedef itk::ImageFileWriter < itk::Image <double, 3> > itkOutputWriter;
	
    typedef itk::ImageFileReader < itk::Image <unsigned char, 3> > itkMaskReader;
    itkMaskReader::Pointer maskRead = itkMaskReader::New();
    maskRead->SetFileName(maskName.c_str());
    maskRead->Update();
    
    typedef itk::VectorImage<double,3> LogTensorImageType;
    typedef itk::ImageFileReader <LogTensorImageType> itkInputReader;
    
    typedef anima::LocalPatchCovarianceDistanceImageFilter<double> MainFilterType;
    
    MainFilterType::Pointer mainFilter = MainFilterType::New();
    mainFilter->SetComputationMask(maskRead->GetOutput());
    mainFilter->SetNumberOfThreads(nbpArg.getValue());

	mainFilter->SetPatchHalfSize(patchHSArg.getValue());
	
    ifstream fileIn(dataLTName.c_str());
    itkInputReader::Pointer ltReader;
    
    unsigned int count = 0;
    while (!fileIn.eof())
    {
        char tmpStr[2048];
        fileIn.getline(tmpStr,2048);
        
        if (strcmp(tmpStr,"") == 0)
            continue;
        
        std::cout << "Loading database image " << tmpStr << "..." << std::endl;
        ltReader = itkInputReader::New();
        ltReader->SetFileName(tmpStr);
        ltReader->Update();
        
        mainFilter->SetInput(count,ltReader->GetOutput());
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
        return 1;
    }
    
    itkOutputWriter::Pointer resultWriter = itkOutputWriter::New();
    resultWriter->SetFileName(resArg.getValue());
	resultWriter->SetUseCompression(true);
	resultWriter->SetInput(mainFilter->GetOutput(0));
    
    resultWriter->Update();
    
    if (resStdArg.getValue() != "")
    {
        itkOutputWriter::Pointer resultStdWriter = itkOutputWriter::New();
        resultStdWriter->SetFileName(resStdArg.getValue());
        resultStdWriter->SetUseCompression(true);
        resultStdWriter->SetInput(mainFilter->GetOutput(1));
        
        resultStdWriter->Update();
    }
    
    return 0;
}
