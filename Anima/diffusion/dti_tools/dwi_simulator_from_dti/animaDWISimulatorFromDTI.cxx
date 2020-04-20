#include <animaDWISimulatorFromDTIImageFilter.h>

#include <itkTimeProbe.h>

#include <animaGradientFileReader.h>
#include <animaReadWriteFunctions.h>

#include <tclap/CmdLine.h>

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
    
    FilterType::Pointer mainFilter = FilterType::New();
    
    typedef anima::GradientFileReader < std::vector <float>, float > GFReaderType;
    GFReaderType gfReader;
    gfReader.SetGradientFileName(gradsArg.getValue());
    gfReader.SetBValueBaseString(bvalArg.getValue());
    gfReader.SetGradientIndependentNormalization(false);
    gfReader.Update();

    std::vector< std::vector<float> > directions = gfReader.GetGradients();
    
    for(unsigned int i = 0;i < directions.size();++i)
        mainFilter->AddGradientDirection(i,directions[i]);
    
    std::vector <float> mb = gfReader.GetBValues();
    mainFilter->SetBValuesList(mb);

    typedef itk::ImageFileReader <FilterType::S0ImageType> S0ImageReaderType;
    FilterType::S0ImageType::Pointer s0Image;
    try
    {
        s0Image = anima::readImage <FilterType::S0ImageType> (s0ValueArg.getValue());
        mainFilter->SetS0Image(s0Image);
    }
    catch(itk::ExceptionObject &e)
    {
        mainFilter->SetS0Value(std::stod(s0ValueArg.getValue()));
    }

    mainFilter->SetInput(anima::readImage <InputImageType> (inArg.getValue()));
    mainFilter->SetNumberOfWorkUnits(nbpArg.getValue());
    
    itk::TimeProbe tmpTimer;
    tmpTimer.Start();
    
    try
    {
        mainFilter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }

    tmpTimer.Stop();
    
    std::cout << "Simulation done in " << tmpTimer.GetTotal() << " s" << std::endl;
    std::cout << "Writing result to : " << resArg.getValue() << std::endl;

    anima::writeImage <Image4DType> (resArg.getValue(), mainFilter->GetOutputAs4DImage());

    return EXIT_SUCCESS;
}
