#include <sstream>
#include <iomanip>

#include <tclap/CmdLine.h>

#include <itksys/SystemTools.hxx>

#include <animaNoiseGeneratorImageFilter.h>
#include <animaReadWriteFunctions.h>
#include <animaRetrieveImageTypeMacros.h>
#include <animaKummerFunctions.h>

struct arguments
{
    std::string input, output;
    double sigma;
    unsigned int nreplicates, nthreads;
    bool gaussianNoise;
};

template <class ComponentType, unsigned int InputDim>
void
addNoise(itk::ImageIOBase::Pointer imageIO, const arguments &args)
{
    typedef itk::Image<ComponentType,InputDim> InputImageType;
    typedef anima::NoiseGeneratorImageFilter<InputImageType> NoiseFilterType;
    typedef typename NoiseFilterType::OutputImageType OutputImageType;
    
    typename NoiseFilterType::Pointer mainFilter = NoiseFilterType::New();
    mainFilter->SetNumberOfWorkUnits(args.nthreads);
    mainFilter->SetInput(anima::readImage<InputImageType>(args.input));
    mainFilter->SetNoiseSigma(args.sigma);
    mainFilter->SetUseGaussianDistribution(args.gaussianNoise);
    mainFilter->SetNumberOfReplicates(args.nreplicates);
    mainFilter->Update();
    
    if (args.nreplicates == 1)
        anima::writeImage<OutputImageType>(args.output, mainFilter->GetOutput(0));
    else
    {
        std::string filePrefix = itksys::SystemTools::GetFilenameWithoutExtension(args.output);
        std::string fileExtension = itksys::SystemTools::GetFilenameExtension(args.output);
        
        for (unsigned int i = 0;i < args.nreplicates;++i)
        {
            std::stringstream ss;
            ss << std::setw(3) << std::setfill('0') << (i + 1);
            std::string tmpStr = ss.str();
            anima::writeImage<OutputImageType>(filePrefix + tmpStr + fileExtension, mainFilter->GetOutput(i));
        }
    }
}

template <class ComponentType >
void
retrieveNbDimensions(itk::ImageIOBase::Pointer imageIO, const arguments &args)
{
    ANIMA_RETRIEVE_NUMBER_OF_DIMENSIONS(imageIO, ComponentType, addNoise, imageIO, args)
}

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","inputfile","Input image on which to add noise.",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","Output image. If number of replicates is 1, it is used as provided for saving the noisy dataset with the corresponding file name. If number of replicates > 1, the extension is extracted as desired output file format and the file name without extension is used as base file name for all replicates.",true,"","output file prefix",cmd);
    TCLAP::ValueArg<std::string> refArg("r","reference","Masked reference image over which the average SNR is defined.",true,"","reference image",cmd);
    
    TCLAP::ValueArg<double> snrArg("s","snr","Average SNR in dB over the reference image (default: 25dB).",false,25.0,"mean snr",cmd);
    TCLAP::ValueArg<unsigned int> repArg("n","num-replicates","Number of independent noisy datasets to generate (default: 1).",false,1,"number of replicates",cmd);
    
    TCLAP::SwitchArg gaussArg("G","gauss-noise","Adds Gaussian noise instead of Rician noise.",cmd,false);
    TCLAP::SwitchArg matchArg("M","match-snr","Make Gaussian noise comparable to Rician noise in terms of SNR.",cmd,false);
    TCLAP::SwitchArg verboseArg("V","verbose","Outputs noise calculations to the console.",cmd,false);

    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default : all cores).",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

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
    
    // Find average signal value in reference image
    typedef itk::Image<float,3> Input3DImageType;
    typedef itk::ImageRegionConstIterator<Input3DImageType> Input3DIteratorType;
    
    Input3DImageType::Pointer referenceImage = anima::readImage<Input3DImageType>(refArg.getValue());
    Input3DIteratorType refItr(referenceImage, referenceImage->GetLargestPossibleRegion());
    
    double meanReferenceValue = 0;
    unsigned int numNonZeroVoxels = 0;
    
    while (!refItr.IsAtEnd())
    {
        double referenceValue = refItr.Get();
        
        if (referenceValue == 0)
        {
            ++refItr;
            continue;
        }
        
        meanReferenceValue *= (numNonZeroVoxels / (numNonZeroVoxels + 1.0));
        meanReferenceValue += referenceValue / (numNonZeroVoxels + 1.0);
        ++numNonZeroVoxels;
        ++refItr;
    }
    
    // Retrieve Rice sigma parameter
    double riceSigmaParameter = meanReferenceValue / snr;
    
    // Compute mean and variance from Rice distribution
    double snrValue = meanReferenceValue / riceSigmaParameter;
    double laguerreArgument = -1.0 * snrValue * snrValue / 2.0;
    double laguerreValue = anima::KummerFunction(laguerreArgument, -0.5, 1.0);
    double riceMeanValue = riceSigmaParameter * std::sqrt(M_PI / 2.0) * laguerreValue;
    double riceStdValue = 2.0 * riceSigmaParameter * riceSigmaParameter + meanReferenceValue * meanReferenceValue - M_PI * riceSigmaParameter * riceSigmaParameter * laguerreValue * laguerreValue / 2.0;
    if (riceStdValue < 0.0)
    {
        std::cerr << "The evaluation of the Kummer M function was too imprecise and produced a negative variance." << std::endl;
        return EXIT_FAILURE;
    }
    riceStdValue = std::sqrt(riceStdValue);
    
    // Deduce effective SNR
    double effSnrValue = riceMeanValue / riceStdValue;
    
    // Compute Gaussian sigma parameter
    double gaussSigmaParameter = (matchArg.isSet()) ? meanReferenceValue / effSnrValue : riceSigmaParameter;
    
    // Set noise standard deviation
    double noiseSigmaParameter = (gaussArg.isSet()) ? gaussSigmaParameter : riceSigmaParameter;
    
    if (verboseArg.isSet())
    {
        // Output some information about the noise structure in the console
        std::cout << " - User-defined average SNR: " << snrArg.getValue() << " dB." << std::endl;
        std::cout << " - Reference image used for noise variance calculation: " << refArg.getValue() << std:: endl;
        std::cout << " - Average signal in foreground voxels of the reference image: " << meanReferenceValue << std::endl;
        std::cout << " - Resulting coil standard deviation (Rice sigma parameter): " << riceSigmaParameter << std::endl;
        std::cout << " - Mean value of Rice distribution: " << riceMeanValue << std::endl;
        std::cout << " - Standard deviation of Rice distribution: " << riceStdValue << std::endl;
        std::cout << " - Effective average SNR: " << 20.0 * std::log10(effSnrValue) << " dB." << std::endl;
        
        if (gaussArg.isSet())
            std::cout << " - Noise distribution: Gaussian" << std::endl;
        else
            std::cout << " - Noise distribution: Rician" << std::endl;
        
        std::cout << " - Resulting noise sigma parameter: " << noiseSigmaParameter << std::endl;
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
    
    arguments args;
    
    args.input = inArg.getValue();
    args.output = outArg.getValue();
    args.sigma = noiseSigmaParameter;
    args.nreplicates = repArg.getValue();
    args.nthreads = nbpArg.getValue();
    args.gaussianNoise = gaussArg.isSet();
    
    try
    {
        ANIMA_RETRIEVE_COMPONENT_TYPE(imageIO, retrieveNbDimensions, imageIO, args);
    }
    catch (itk::ExceptionObject &err)
    {
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}
