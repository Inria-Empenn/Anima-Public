#include <animaBootstrap4DVolumeImageFilter.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkTimeProbe.h>
#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","inputlist","File containing a list of 4D images",true,"","input 4D images",cmd);
    TCLAP::ValueArg<std::string> infoArg("l","infolist","File containing a list of info on each 4D image (one info line per sub-3D volume)",false,"","info list on 4D images",cmd);
    TCLAP::ValueArg<std::string> resArg("o","outputfile","Result sampled 4D image",true,"","result 4D image",cmd);
    TCLAP::ValueArg<std::string> resInfoArg("O","outputinfo","Result information gathered from input information",false,"","result info text file",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    typedef anima::Bootstrap4DVolumeImageFilter <double> MainFilterType;

    typedef itk::ImageFileReader < MainFilterType::TInputImage > ImageReaderType;
    typedef itk::ImageFileWriter < MainFilterType::TOutputImage > ImageWriterType;

    MainFilterType::Pointer mainFilter = MainFilterType::New();

    std::ifstream inputFile(inArg.getValue().c_str());

    if (!inputFile.is_open())
    {
        std::cerr << "Please provide usable file with input files" << std::endl;
        return -1;
    }

    unsigned int numInput = 0;
    while (!inputFile.eof())
    {
        char tmpStr[2048];
        inputFile.getline(tmpStr,2048);

        if (strcmp(tmpStr,"") == 0)
            continue;

        ImageReaderType::Pointer tmpReader = ImageReaderType::New();
        tmpReader->SetFileName(tmpStr);
        tmpReader->Update();
        mainFilter->SetInput(numInput,tmpReader->GetOutput());
        numInput++;
    }

    if (infoArg.getValue() != "")
    {
        std::ifstream infoList(infoArg.getValue().c_str());

        unsigned int i = 0;
        while (!infoList.eof())
        {
            char tmpStr[2048];
            infoList.getline(tmpStr,2048);

            if (strcmp(tmpStr,"") == 0)
                continue;

            mainFilter->SetNthInformationFile(i,tmpStr);
            ++i;
        }

        infoList.close();
    }

    mainFilter->Update();

    std::cout << "Writing result to : " << resArg.getValue() << std::endl;

    ImageWriterType::Pointer outWriter = ImageWriterType::New();
    outWriter->SetInput(mainFilter->GetOutput());
    outWriter->SetFileName(resArg.getValue());
    outWriter->SetUseCompression(true);

    outWriter->Update();

    if ((resInfoArg.getValue() != "")&&(mainFilter->GetOutputInformation().size() != 0))
    {
        std::ofstream outInfo(resInfoArg.getValue().c_str());
        for (unsigned int i = 0;i < mainFilter->GetOutputInformation().size();++i)
            outInfo << mainFilter->GetOutputInformation()[i] << std::endl;

        outInfo.close();
    }

    return 0;
}
