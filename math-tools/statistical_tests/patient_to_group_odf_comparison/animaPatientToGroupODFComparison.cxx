#include <iostream>
#include <tclap/CmdLine.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <animaPatientToGroupODFComparisonImageFilter.h>
#include <itkTimeProbe.h>

#include <animaODFFunctions.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',"1.0");
	
    TCLAP::ValueArg<std::string> refODFArg("i","inputodf","ODF Test Image",true,"","ODF test image",cmd);
	
    TCLAP::ValueArg<std::string> dataODFArg("I","databaseodf","ODF Database Image List",true,"","ODF database image list",cmd);
    
	TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask",true,"","computation mask",cmd);
    TCLAP::ValueArg<std::string> resArg("o","outputname","Z-Score output image",true,"","Z-Score output image",cmd);
    TCLAP::ValueArg<std::string> resPValArg("O","outpvalname","P-value output image",true,"","P-Value output image",cmd);

    TCLAP::ValueArg<std::string> statTestArg("t","stat-test","Statistical test to use ([fisher],chi)",false,"fisher","statistical test",cmd);
    TCLAP::ValueArg<double> expVarArg("e","expvar","PCA threshold: threshold on eigenvalues to compute the new basis (default: 0.5)",false,0.5,"PCA threshold",cmd);
	TCLAP::ValueArg<unsigned int> numEigenArg("E","numeigenpca","Number of eigenvalues to keep (default: 6)",false,6,"Number of PCA eigen values",cmd);
    
    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default: all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
	
	TCLAP::ValueArg<unsigned int> nbThetaArg("T","theta","Number of theta values (theta varies between 0 and pi/2",false,0,"number of theta values",cmd);
    TCLAP::ValueArg<unsigned int> nbPhiArg("P","phi","Number of phi values (theta varies between 0 and 2 pi",false,0,"number of phi values",cmd);
	
	TCLAP::ValueArg<std::string> samplesFileNameArg("d","sampledirectionsfile","Samples directions in a text file",false,"","Samples directions file",cmd);
	
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
	
    typedef itk::ImageFileWriter < itk::Image <double, 3> > itkOutputWriter;
	
    typedef itk::ImageFileReader < itk::Image <unsigned char, 3> > itkMaskReader;
    itkMaskReader::Pointer maskRead = itkMaskReader::New();
    maskRead->SetFileName(maskName.c_str());
    maskRead->Update();
	
    typedef itk::VectorImage<double,3> ODFImageType;
    typedef itk::ImageFileReader <ODFImageType> itkInputReader;
	
    typedef anima::PatientToGroupODFComparisonImageFilter<double> MAOZScoreImageFilterType;
    
    MAOZScoreImageFilterType::Pointer mainFilter = MAOZScoreImageFilterType::New();
    mainFilter->SetComputationMask(maskRead->GetOutput());
    mainFilter->SetNumberOfThreads(nbpArg.getValue());
	
	itkInputReader::Pointer testODFReader = itkInputReader::New();
	testODFReader->SetFileName(refODFArg.getValue());
	testODFReader->Update();
	mainFilter->SetInput(testODFReader->GetOutput());
	
	mainFilter->SetExplainedRatio(expVarArg.getValue());
	mainFilter->SetNumEigenValuesPCA(numEigenArg.getValue());
	
    std::vector < std::vector <double> > sampleDirections = anima::InitializeSampleDirections(nbThetaArg.getValue(),nbPhiArg.getValue(),
                                                                                              samplesFileNameArg.getValue());
	
	mainFilter->SetSampleDirections(sampleDirections);
	
    if (statTestArg.getValue() == "chi")
        mainFilter->SetStatisticalTestType(MAOZScoreImageFilterType::CHI_SQUARE);
    else
        mainFilter->SetStatisticalTestType(MAOZScoreImageFilterType::FISHER);

    std::ifstream fileIn(dataODFName.c_str());
    if (!fileIn.is_open())
    {
        std::cerr << "Could not open ODF data file (" << dataODFName << ")" << std::endl;
        return EXIT_FAILURE;
    }

    itkInputReader::Pointer odfReader;
    
    while (!fileIn.eof())
    {
        char tmpStr[2048];
        fileIn.getline(tmpStr,2048);
        
        if (strcmp(tmpStr,"") == 0)
            continue;
        
        std::cout << "Loading ODF image " << tmpStr << "..." << std::endl;
        odfReader = itkInputReader::New();
        odfReader->SetFileName(tmpStr);
        odfReader->Update();
        
        mainFilter->AddDatabaseInput(odfReader->GetOutput());
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
	
    itkOutputWriter::Pointer resultPValWriter = itkOutputWriter::New();
    resultPValWriter->SetFileName(resPValArg.getValue());
    resultPValWriter->SetUseCompression(true);
    resultPValWriter->SetInput(mainFilter->GetOutput(1));

    resultPValWriter->Update();

    return 0;
}
