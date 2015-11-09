#include <iostream>
#include <tclap/CmdLine.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <animaPatientToGroupComparisonImageFilter.h>
#include <itkTimeProbe.h>

using namespace std;

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> refLTArg("i","input","Test Image",true,"","test image",cmd);
    TCLAP::ValueArg<std::string> dataLTArg("I","database","Database Image List",true,"","database image list",cmd);
    
	TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask",true,"","computation mask",cmd);
    TCLAP::ValueArg<std::string> resArg("o","outputname","Z-Score output image",true,"","Z-Score output image",cmd);
    TCLAP::ValueArg<std::string> resPValArg("O","outpvalname","P-value output image",true,"","P-Value output image",cmd);

    TCLAP::ValueArg<std::string> statTestArg("t","stat-test","Statistical test to use ([fisher],chi)",false,"fisher","statistical test",cmd);
    TCLAP::ValueArg<double> expVarArg("e","expvar","PCA threshold: threshold on eigenvalues to compute the new basis (default: 0.5)",false,0.5,"PCA threshold",cmd);
    TCLAP::ValueArg<unsigned int> numEigenArg("E","numeigenpca","Number of eigenvalues to keep (default: 6)",false,6,"Number of PCA eigen values",cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default: all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
	
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
    
    typedef anima::PatientToGroupComparisonImageFilter<double> MAZScoreImageFilterType;
    
    MAZScoreImageFilterType::Pointer mainFilter = MAZScoreImageFilterType::New();
    mainFilter->SetComputationMask(maskRead->GetOutput());
    mainFilter->SetNumberOfThreads(nbpArg.getValue());
	
	itkInputReader::Pointer testLTReader = itkInputReader::New();
	testLTReader->SetFileName(refLTArg.getValue());
	testLTReader->Update();
	mainFilter->SetInput(testLTReader->GetOutput());

    mainFilter->SetExplainedRatio(expVarArg.getValue());
    mainFilter->SetNumEigenValuesPCA(numEigenArg.getValue());

    if (statTestArg.getValue() == "chi")
        mainFilter->SetStatisticalTestType(MAZScoreImageFilterType::CHI_SQUARE);
    else
        mainFilter->SetStatisticalTestType(MAZScoreImageFilterType::FISHER);

    ifstream fileIn(dataLTName.c_str());
    if (!fileIn.is_open())
    {
        std::cerr << "Could not open ODF data file (" << dataLTName << ")" << std::endl;
        return EXIT_FAILURE;
    }
    
    itkInputReader::Pointer ltReader;
    
    while (!fileIn.eof())
    {
        char tmpStr[2048];
        fileIn.getline(tmpStr,2048);
        
        if (strcmp(tmpStr,"") == 0)
            continue;
        
        std::cout << "Loading tensor image " << tmpStr << "..." << std::endl;
        ltReader = itkInputReader::New();
        ltReader->SetFileName(tmpStr);
        ltReader->Update();
        
        mainFilter->AddDatabaseInput(ltReader->GetOutput());
    }
    fileIn.close();
  		
    try
    {
        mainFilter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return 1;
    }
    
    itkOutputWriter::Pointer resultWriter = itkOutputWriter::New();
    resultWriter->SetFileName(resArg.getValue());
    resultWriter->SetUseCompression(true);
    resultWriter->SetInput(mainFilter->GetOutput(0));

    resultWriter->Update();

    itkOutputWriter::Pointer resultPValWriter = itkOutputWriter::New();
    resultPValWriter->SetFileName(resPValArg.getValue());
    resultPValWriter->SetUseCompression(true);
    resultPValWriter->SetInput(mainFilter->GetOutput(1));

    resultPValWriter->Update();

    return 0;
}
