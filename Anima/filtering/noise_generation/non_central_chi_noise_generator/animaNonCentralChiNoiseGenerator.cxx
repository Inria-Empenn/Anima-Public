#include <tclap/CmdLine.h>

#include <animaNonCentralChiNoiseGeneratorImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

using namespace std;

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',"1.0");

    TCLAP::ValueArg<std::string> inArg("i","inputfile","Input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","Output image",true,"","output image",cmd);

    TCLAP::ValueArg<double> backVarArg("b","backgroundvariance","Background variance (default: 20)",false,20,"background variance",cmd);
    TCLAP::ValueArg<double> relStdArg("s","std","Desired inverse SNR (default: 0.05)",false,0.05,"noise std",cmd);

    TCLAP::ValueArg<double> refValArg("r","ref-val","Reference value on which relative noise is applied (default: computed from image)",false,0,"reference value",cmd);
    TCLAP::SwitchArg avg4dMeanValArg("A","average-4d-mean","Automatically computed mean value of tissue is over all 4D volume",cmd,false);

    TCLAP::ValueArg<unsigned int> numCoilsArg("c","num-coils","Number of coils (default : 1 = Rician noise)",false,1,"number of coils",cmd);
    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default : all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    string refName, resName;
    refName = inArg.getValue();
    resName = outArg.getValue();

    double backVar = backVarArg.getValue();
    double relStd = relStdArg.getValue();

    // Find out the type of the image in file
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inArg.getValue().c_str(),itk::ImageIOFactory::ReadMode);

    if( !imageIO )
    {
        std::cerr << "Itk could not find a suitable IO factory for the input" << std::endl;
        return EXIT_FAILURE;
    }

    // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
    imageIO->SetFileName(inArg.getValue());
    imageIO->ReadImageInformation();

    const unsigned int nbDimension = imageIO->GetNumberOfDimensions ();

    switch (nbDimension)
    {
        case 4:
        {
            typedef itk::Image<float,4> ImageType;
            typedef itk::ImageFileReader <ImageType> itkInputReader;
            typedef itk::ImageFileWriter <ImageType> itkOutputWriter;

            typedef anima::NonCentralChiNoiseGeneratorImageFilter<4> NoiseFilterType;

            NoiseFilterType::Pointer mainFilter = NoiseFilterType::New();
            mainFilter->SetNumberOfThreads(nbpArg.getValue());

            itkInputReader::Pointer imageReader = itkInputReader::New();
            imageReader->SetFileName(refName.c_str());
            imageReader->Update();

            mainFilter->SetInput(imageReader->GetOutput());
            mainFilter->SetBackgroundVariance(backVar);
            mainFilter->SetRefVal(refValArg.getValue());
            mainFilter->SetRelStdDevGaussianNoise(relStd);
            mainFilter->SetNumberOfCoils(numCoilsArg.getValue());

            // If a reference value is specified, ignore average flag
            if (refValArg.getValue() == 0)
                mainFilter->SetAverageMeanValOnAllVolumes(avg4dMeanValArg.isSet());
            else
                mainFilter->SetAverageMeanValOnAllVolumes(true);

            mainFilter->Update();

            itkOutputWriter::Pointer resultWriter = itkOutputWriter::New();
            resultWriter->SetFileName(resName.c_str());
            resultWriter->SetUseCompression(true);
            resultWriter->SetInput(mainFilter->GetOutput());

            resultWriter->Update();

            break;
        }

        case 3:
        default:
        {
            typedef itk::Image<float,3> ImageType;
            typedef itk::ImageFileReader <ImageType> itkInputReader;
            typedef itk::ImageFileWriter <ImageType> itkOutputWriter;

            typedef anima::NonCentralChiNoiseGeneratorImageFilter<3> NoiseFilterType;

            NoiseFilterType::Pointer mainFilter = NoiseFilterType::New();
            mainFilter->SetNumberOfThreads(nbpArg.getValue());

            itkInputReader::Pointer imageReader = itkInputReader::New();
            imageReader->SetFileName(refName.c_str());
            imageReader->Update();

            mainFilter->SetInput(imageReader->GetOutput());
            mainFilter->SetBackgroundVariance(backVar);
            mainFilter->SetRefVal(refValArg.getValue());
            mainFilter->SetRelStdDevGaussianNoise(relStd);
            mainFilter->SetNumberOfCoils(numCoilsArg.getValue());
            mainFilter->SetAverageMeanValOnAllVolumes(true);

            mainFilter->Update();

            itkOutputWriter::Pointer resultWriter = itkOutputWriter::New();
            resultWriter->SetFileName(resName.c_str());
            resultWriter->SetUseCompression(true);
            resultWriter->SetInput(mainFilter->GetOutput());

            resultWriter->Update();

            break;
        }
    }

    return 0;
}
