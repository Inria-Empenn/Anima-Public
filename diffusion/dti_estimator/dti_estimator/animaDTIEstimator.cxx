#include <animaDTIEstimationImageFilter.h>
#include <animaGradientFileReader.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <animaReadWriteFunctions.h>
#include <itkTimeProbe.h>

#include <tclap/CmdLine.h>

int main(int argc,  char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',"1.0");

    TCLAP::ValueArg<std::string> inArg("i","inputdwi","DWI volume",true,"","DWI volume",cmd);
    TCLAP::ValueArg<std::string> resArg("o","output","Result DTI image",true,"","result DTI image",cmd);
    TCLAP::ValueArg<std::string> b0OutArg("O","output-b0","Resulting B0 image",true,"","result B0 image",cmd);

    TCLAP::ValueArg<std::string> gradsArg("g","grad","Input gradients",true,"","Input gradients",cmd);
    TCLAP::ValueArg<std::string> bvalArg("b","bval","Input b-values",true,"","Input b-values",cmd);

    TCLAP::SwitchArg keepDegArg("K","keep-degenerated","Keep degenerated values",cmd,false);

    TCLAP::ValueArg<unsigned int> b0ThrArg("t","b0thr","B0 threshold (default : 0)",false,0,"B0 threshold for computation",cmd);
    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default: all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try{
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    typedef anima::DTIEstimationImageFilter <float> FilterType;
    typedef FilterType::InputImageType InputImageType;
    typedef FilterType::Image4DType Image4DType;
    typedef FilterType::OutputImageType VectorImageType;

    typedef itk::ImageFileWriter <InputImageType> OutputB0ImageWriterType;
    typedef itk::ImageFileWriter <VectorImageType> VectorImageWriterType;

    FilterType::Pointer mainFilter = FilterType::New();

    typedef anima::GradientFileReader < std::vector < float >, float > GFReaderType;
    GFReaderType gfReader;
    gfReader.SetGradientFileName(gradsArg.getValue());
    gfReader.SetBValueBaseString(bvalArg.getValue());
    gfReader.SetGradientIndependentNormalization(false);

    gfReader.Update();

    GFReaderType::GradientVectorType directions = gfReader.GetGradients();

    for(unsigned int i = 0;i < directions.size();++i)
        mainFilter->AddGradientDirection(i,directions[i]);

    GFReaderType::BValueVectorType mb = gfReader.GetBValues();

    mainFilter->SetBValuesList(mb);

    anima::setMultipleImageFilterInputsFromFileName<InputImageType,Image4DType,FilterType>(inArg.getValue(), mainFilter);

    mainFilter->SetB0Threshold(b0ThrArg.getValue());
    mainFilter->SetNumberOfThreads(nbpArg.getValue());
    mainFilter->SetRemoveDegeneratedTensors(!keepDegArg.isSet());

    itk::TimeProbe tmpTimer;

    tmpTimer.Start();

    try
    {
        mainFilter->Update();
    } catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return -1;
    }

    tmpTimer.Stop();

    std::cout << "Estimation done in " << tmpTimer.GetTotal() << " s" << std::endl;

    std::cout << "Writing result to : " << resArg.getValue() << std::endl;

    VectorImageType::Pointer output = mainFilter->GetOutput();
    output->DisconnectPipeline();

    VectorImageWriterType::Pointer vectorWriter = VectorImageWriterType::New();
    vectorWriter->SetInput(output);
    vectorWriter->SetFileName(resArg.getValue());
    vectorWriter->SetUseCompression(true);

    vectorWriter->Update();

    if (b0OutArg.getValue() != "")
    {
        OutputB0ImageWriterType::Pointer b0Writer = OutputB0ImageWriterType::New();
        b0Writer->SetInput(mainFilter->GetEstimatedB0Image());
        b0Writer->SetFileName(b0OutArg.getValue());
        b0Writer->SetUseCompression(true);

        b0Writer->Update();
    }

    return 0;
}
