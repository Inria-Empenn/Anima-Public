#include <sstream>
#include <iomanip>

#include <tclap/CmdLine.h>

#include <itksys/SystemTools.hxx>

#include <animaNoiseGeneratorImageFilter.h>
#include <animaReadWriteFunctions.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","inputfile","Input image on which to add noise.",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","Output image. If number of replicates is 1, it is used as provided for saving the noisy dataset with the corresponding file name. If number of replicates > 1, the extension is extracted as desired output file format and the file name without extension is used as base file name for all replicates.",true,"","output file prefix",cmd);
    TCLAP::ValueArg<std::string> refArg("r","reference","Masked reference image over which the average SNR is defined.",true,"","reference image",cmd);
    
    TCLAP::ValueArg<double> snrArg("s","snr","Average SNR in dB over the reference image (default: 25dB).",false,25.0,"mean snr",cmd);
    TCLAP::ValueArg<unsigned int> repArg("n","num-replicates","Number of independent noisy datasets to generate (default: 1).",false,1,"number of replicates",cmd);
    
    TCLAP::SwitchArg gaussArg("G","gauss-noise","Adds Gaussian noise instead of Rician noise.",cmd,false);
    TCLAP::SwitchArg verboseArg("V","verbose","Outputs noise calculations to the console.",cmd,false);

    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default : all cores).",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    // Put SNR back in linear scale
    double snr = std::pow(10.0, snrArg.getValue() / 20.0);
    
    // Find average B0 value
    typedef itk::Image<float,3> Input3DImageType;
    typedef itk::ImageRegionConstIterator<Input3DImageType> Iterator3DImageType;
    
    Input3DImageType::Pointer b0Image = anima::readImage<Input3DImageType>(refArg.getValue());
    Iterator3DImageType b0Itr(b0Image, b0Image->GetLargestPossibleRegion());
    
    double meanB0Value = 0;
    unsigned int numNonZeroVoxels = 0;
    
    while (!b0Itr.IsAtEnd())
    {
        double b0Value = b0Itr.Get();
        
        if (b0Value == 0)
        {
            ++b0Itr;
            continue;
        }
        
        meanB0Value *= (numNonZeroVoxels / (numNonZeroVoxels + 1.0));
        meanB0Value += b0Value / (numNonZeroVoxels + 1.0);
        ++numNonZeroVoxels;
        ++b0Itr;
    }
    
    // Retrieve standard deviation of noise
    double stdDev = meanB0Value / snr;
    
    if (verboseArg.isSet())
    {
        // Output some information about the noise structure in the console
        std::cout << " - User-defined average SNR: " << snrArg.getValue() << " dB." << std::endl;
        std::cout << " - Average signal in foreground voxels of the input B0 image: " << meanB0Value << std::endl;
        if (gaussArg.isSet())
            std::cout << " - Noise distribution: Gaussian" << std::endl;
        else
            std::cout << " - Noise distribution: Rician" << std::endl;
        std::cout << " - Resulting standard deviation for noise simulation: " << stdDev << std::endl;
    }

    // Find out the type of the image in file
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inArg.getValue().c_str(),itk::ImageIOFactory::ReadMode);

    if (!imageIO)
    {
        std::cerr << "Itk could not find a suitable IO factory for the input" << std::endl;
        return EXIT_FAILURE;
    }

    // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
    imageIO->SetFileName(inArg.getValue());
    imageIO->ReadImageInformation();

    const unsigned int nbDimension = imageIO->GetNumberOfDimensions();

    switch (nbDimension)
    {
        case 4:
        {
            typedef itk::Image<float,4> ImageType;
            typedef anima::NoiseGeneratorImageFilter<4> NoiseFilterType;

            NoiseFilterType::Pointer mainFilter = NoiseFilterType::New();
            mainFilter->SetNumberOfThreads(nbpArg.getValue());
            mainFilter->SetInput(anima::readImage<ImageType>(inArg.getValue()));
            mainFilter->SetStandardDeviation(stdDev);
            mainFilter->SetUseGaussianDistribution(gaussArg.isSet());
            mainFilter->SetNumberOfReplicates(repArg.getValue());
            mainFilter->Update();

            if (repArg.getValue() == 1)
                anima::writeImage<ImageType>(outArg.getValue(), mainFilter->GetOutput(0));
            else
            {
                std::string filePrefix = itksys::SystemTools::GetFilenameWithoutExtension(outArg.getValue());
                std::string fileExtension = itksys::SystemTools::GetFilenameExtension(outArg.getValue());
                
                for (unsigned int i = 0;i < repArg.getValue();++i)
                {
                    std::stringstream ss;
                    ss << std::setw(3) << std::setfill('0') << (i + 1);
                    std::string tmpStr = ss.str();
                    anima::writeImage<ImageType>(filePrefix + tmpStr + fileExtension, mainFilter->GetOutput(i));
                }
            }
            
            break;
        }

        case 3:
        default:
        {
            typedef itk::Image<float,3> ImageType;
            typedef anima::NoiseGeneratorImageFilter<3> NoiseFilterType;

            NoiseFilterType::Pointer mainFilter = NoiseFilterType::New();
            mainFilter->SetNumberOfThreads(nbpArg.getValue());
            mainFilter->SetInput(anima::readImage<ImageType>(inArg.getValue()));
            mainFilter->SetStandardDeviation(stdDev);
            mainFilter->SetUseGaussianDistribution(gaussArg.isSet());
            mainFilter->SetNumberOfReplicates(repArg.getValue());
            mainFilter->Update();

            if (repArg.getValue() == 1)
                anima::writeImage<ImageType>(outArg.getValue(), mainFilter->GetOutput(0));
            else
            {
                std::string filePrefix = itksys::SystemTools::GetFilenameWithoutExtension(outArg.getValue());
                std::string fileExtension = itksys::SystemTools::GetFilenameExtension(outArg.getValue());
                
                for (unsigned int i = 0;i < repArg.getValue();++i)
                {
                    std::stringstream ss;
                    ss << std::setw(3) << std::setfill('0') << (i + 1);
                    std::string tmpStr = ss.str();
                    anima::writeImage<ImageType>(filePrefix + tmpStr + fileExtension, mainFilter->GetOutput(i));
                }
            }
            
            break;
        }
    }

    return EXIT_SUCCESS;
}
