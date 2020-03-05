#include <animaDWISimulatorFromDTIImageFilter.h>

#include <fstream>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbe.h>

#include <tclap/CmdLine.h>

void loadGradients(std::string &gradsArg, std::vector< std::vector<float> >& directions)
{
    std::ifstream gradFile(gradsArg.c_str());
    
    if (!gradFile.is_open())
    {
        std::cerr << "Please provide file with gradient information" << std::endl;
        abort();
    }
    
    std::vector <float> gradTmp(3,0);
    directions.clear();
    
    while (!gradFile.eof())
    {
        char tmpStr[2048];
        gradFile.getline(tmpStr,2048);
        
        if (strcmp(tmpStr,"") == 0)
            continue;
        
        sscanf(tmpStr,"%f %f %f",&gradTmp[0],&gradTmp[1],&gradTmp[2]);
        
        float normTmp = std::sqrt(gradTmp[0]*gradTmp[0] + gradTmp[1]*gradTmp[1] + gradTmp[2]*gradTmp[2]);
        
        if (normTmp > 0)
        {
            for (unsigned int i = 0;i < 3;++i)
                gradTmp[i] /= normTmp;
        }
        
        directions.push_back(gradTmp);
    }
    
    gradFile.close();
}

void loadBValues(std::vector<std::vector<float> >& directions, std::string &bvalsArg, std::vector<float>& mb)
{
    mb.clear();
    
    std::ifstream bvalFile(bvalsArg.c_str());
    
    if (!bvalFile.is_open())
    {        
        for (unsigned int i = 0;i < directions.size();++i)
            mb.push_back(std::stod(bvalsArg));
        
        return;
    }
    
    float bvalTmp = 0;
    
    while (!bvalFile.eof())
    {
        char tmpStr[2048];
        bvalFile.getline(tmpStr,2048);
        
        if (strcmp(tmpStr,"") == 0)
            continue;
        
        sscanf(tmpStr,"%f",&bvalTmp);
        mb.push_back(bvalTmp);
    }
    
    bvalFile.close();
}

int main(int argc,  char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
	
    TCLAP::ValueArg<std::string> inArg("i","inputdti","DTI volume",true,"","DTI volume",cmd);
	TCLAP::ValueArg<std::string> resArg("o","output","Result DWI image",true,"","result DWI image",cmd);
    TCLAP::ValueArg<std::string> s0ValueArg("s","s0","S0 of DWI (constant value or image)",false,"200","S0 of DWI",cmd);

    TCLAP::ValueArg<std::string> gradsArg("g","grad","Input gradients",true,"","Input gradients",cmd);
    TCLAP::ValueArg<std::string> bvalArg("b","bval","Input b-values",true,"","Input b-values",cmd);

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
    
    typedef anima::DWISimulatorFromDTIImageFilter<float> FilterType;
    typedef FilterType::InputImageType InputImageType;
    typedef FilterType::Image4DType Image4DType;

    typedef itk::ImageFileReader <InputImageType> InputImageReaderType;
    typedef itk::ImageFileWriter <Image4DType> Output4DImageWriterType;
    
    FilterType::Pointer mainFilter = FilterType::New();
    
    std::vector<std::vector<float> > directions;
    loadGradients(gradsArg.getValue(), directions);
    
    for(unsigned int i = 0;i < directions.size();++i)
        mainFilter->AddGradientDirection(i,directions[i]);
    
    std::vector<float> mb;
	loadBValues(directions, bvalArg.getValue(), mb);
    
    mainFilter->SetBValuesList(mb);
    
    typedef itk::ImageFileReader <FilterType::S0ImageType> S0ImageReaderType;
    S0ImageReaderType::Pointer s0TentativeReader = S0ImageReaderType::New();
    s0TentativeReader->SetFileName(s0ValueArg.getValue());

    try
    {
        s0TentativeReader->Update();
        mainFilter->SetS0Image(s0TentativeReader->GetOutput());
    }
    catch(itk::ExceptionObject &e)
    {
        mainFilter->SetS0Value(std::stod(s0ValueArg.getValue()));
    }

    InputImageReaderType::Pointer reader = InputImageReaderType::New();
    reader->SetFileName(inArg.getValue());
    
    reader->Update();
    mainFilter->SetInput(0,reader->GetOutput());

    mainFilter->SetNumberOfWorkUnits(nbpArg.getValue());
    
    itk::TimeProbe tmpTimer;
    
    tmpTimer.Start();
    
    try {
        mainFilter->Update();
    } catch (itk::ExceptionObject &e) {
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }
        
    tmpTimer.Stop();
    
    std::cout << "Simulation done in " << tmpTimer.GetTotal() << " s" << std::endl;

    std::cout << "Writing result to : " << resArg.getValue() << std::endl;

    Output4DImageWriterType::Pointer writer = Output4DImageWriterType::New();
    writer->SetInput(mainFilter->GetOutputAs4DImage());
    writer->SetFileName(resArg.getValue());
    writer->SetUseCompression(true);
    
    writer->Update();

    return EXIT_SUCCESS;
}
