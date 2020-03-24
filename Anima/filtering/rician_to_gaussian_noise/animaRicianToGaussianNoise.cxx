#include <animaRicianToGaussianImageFilter.h>
#include <animaReadWriteFunctions.h>
#include <animaRetrieveImageTypeMacros.h>

#include <tclap/CmdLine.h>

//Update progression of the process
void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
    std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
}

struct arguments
{
    std::string input, output, mask, mean, variance;
    double epsilon, sigma;
    unsigned int maxiters, nthreads;
};

template <class ComponentType>
void
concatenate(const std::vector<typename itk::Image<ComponentType,3>::Pointer> &inputData, const std::string &output)
{
    const unsigned int InputDim = 3;
    typedef itk::Image<ComponentType,InputDim> InputImageType;
    typedef itk::Image<ComponentType,4> OutputImageType;

    typename OutputImageType::Pointer baseImg;
    
    typename OutputImageType::RegionType finalRegion;
    typename OutputImageType::SizeType finalSize;
    typename OutputImageType::IndexType finalIndex;
    typename OutputImageType::PointType finalOrigin;
    typename OutputImageType::SpacingType finalSpacing;
    typename OutputImageType::DirectionType finalDirection;

    for (unsigned int d = 0;d < InputDim;++d)
    {
        finalIndex[d] = inputData[0]->GetLargestPossibleRegion().GetIndex()[d];
        finalSize[d] = inputData[0]->GetLargestPossibleRegion().GetSize()[d];
        finalOrigin[d] = inputData[0]->GetOrigin()[d];
        finalSpacing[d] = inputData[0]->GetSpacing()[d];
        for(unsigned int i = 0;i < InputDim;++i)
            finalDirection[d][i] = inputData[0]->GetDirection()[d][i];
    }
    finalIndex[InputDim] = 0;
    finalSize[InputDim] = 1;
    finalOrigin[InputDim] = 0;
    finalSpacing[InputDim] = 1;

    for(unsigned int i = 0; i < InputDim; ++i)
    {
        finalDirection[InputDim][i] = 0;
        finalDirection[i][InputDim] = 0;
    }
    finalDirection[InputDim][InputDim] = 1;

    finalRegion.SetIndex(finalIndex);
    finalRegion.SetSize(finalSize);

    baseImg = OutputImageType::New();
    baseImg->Initialize();
    baseImg->SetRegions(finalRegion);
    baseImg->SetOrigin(finalOrigin);
    baseImg->SetSpacing(finalSpacing);
    baseImg->SetDirection(finalDirection);
    baseImg->SetNumberOfComponentsPerPixel(1);
    baseImg->Allocate();

    typedef itk::ImageRegionIterator <OutputImageType> FillIteratorType;
    FillIteratorType fillItr(baseImg, baseImg->GetLargestPossibleRegion());
    typedef itk::ImageRegionConstIterator<InputImageType> SourceIteratorType;
    SourceIteratorType srcItr(inputData[0], inputData[0]->GetLargestPossibleRegion());
    while(!srcItr.IsAtEnd())
    {
        fillItr.Set(srcItr.Get());
        ++fillItr; ++srcItr;
    }

    // allocate output
    typename OutputImageType::Pointer outputImg;

    for (unsigned int d = 0; d < InputDim; ++d)
    {
        finalIndex[d] = baseImg->GetLargestPossibleRegion().GetIndex()[d];
        finalSize[d] = baseImg->GetLargestPossibleRegion().GetSize()[d];
        finalOrigin[d] = baseImg->GetOrigin()[d];
        finalSpacing[d] = baseImg->GetSpacing()[d];

        for(unsigned int i = 0;i < InputDim;++i)
            finalDirection[d][i] = baseImg->GetDirection()[d][i];
    }

    for(unsigned int i = 0;i < InputDim + 1;++i)
    {
        finalDirection[InputDim][i] = baseImg->GetDirection()[InputDim][i];
        finalDirection[i][InputDim] = baseImg->GetDirection()[i][InputDim];
    }

    finalIndex[InputDim] = baseImg->GetLargestPossibleRegion().GetIndex()[InputDim];
    finalSize[InputDim] = baseImg->GetLargestPossibleRegion().GetSize()[InputDim] + inputData.size() - 1;
    finalOrigin[InputDim] = baseImg->GetOrigin()[InputDim];
    finalSpacing[InputDim] = baseImg->GetSpacing()[InputDim];

    finalRegion.SetIndex(finalIndex);
    finalRegion.SetSize(finalSize);

    outputImg = OutputImageType::New();
    outputImg->Initialize();
    outputImg->SetRegions(finalRegion);
    outputImg->SetOrigin(finalOrigin);
    outputImg->SetSpacing(finalSpacing);
    outputImg->SetDirection(finalDirection);
    outputImg->SetNumberOfComponentsPerPixel(1);
    outputImg->Allocate();

    // fill with the base:
    typedef itk::ImageRegionIterator <OutputImageType> FillIteratorType;
    fillItr = FillIteratorType(outputImg, outputImg->GetLargestPossibleRegion());
    typedef itk::ImageRegionConstIterator<OutputImageType> BaseIteratorType;
    BaseIteratorType baseItr(baseImg, baseImg->GetLargestPossibleRegion());
    while(!baseItr.IsAtEnd())
    {
        fillItr.Set(baseItr.Get());
        ++fillItr; ++baseItr;
    }

    //fill with the inputs:
    for(unsigned int i = 1; i < inputData.size();++i)
    {
        typedef itk::ImageRegionConstIterator<InputImageType> InputIteratorType;
        InputIteratorType inputItr(inputData[i], inputData[i]->GetLargestPossibleRegion());
        while(!inputItr.IsAtEnd())
        {
            fillItr.Set(inputItr.Get());
            ++fillItr; ++inputItr;
        }
    }

    // write res
    anima::writeImage<OutputImageType>(output, outputImg);
}

template <class ComponentType, unsigned int ImageDimension>
typename itk::Image<ComponentType,ImageDimension>::Pointer
performTransformation(const typename itk::Image<ComponentType,ImageDimension>::Pointer inputImage,
                      const typename itk::Image<ComponentType,ImageDimension>::Pointer meanImage,
                      const typename itk::Image<ComponentType,ImageDimension>::Pointer varImage,
                      const arguments &args)
{
    typedef anima::RicianToGaussianImageFilter<ComponentType,ImageDimension> RiceToGaussianImageFilterType;
    typedef typename RiceToGaussianImageFilterType::InputImageType InputImageType;
    typedef typename RiceToGaussianImageFilterType::OutputImageType OutputImageType;
    typedef typename RiceToGaussianImageFilterType::MaskImageType MaskImageType;
    
    typename MaskImageType::Pointer maskImage;
    if (args.mask != "")
        maskImage = anima::readImage<MaskImageType>(args.mask);
    
    itk::CStyleCommand::Pointer callback;
    typename RiceToGaussianImageFilterType::Pointer mainFilter;
    
    // Estimate scale parameter using SNR fixed point strategy
    mainFilter = RiceToGaussianImageFilterType::New();
    mainFilter->SetInput(inputImage);
    mainFilter->SetMaximumNumberOfIterations(args.maxiters);
    mainFilter->SetEpsilon(args.epsilon);
    mainFilter->SetSigma(args.sigma);
    mainFilter->SetNumberOfWorkUnits(args.nthreads);
    
    if (maskImage)
        mainFilter->SetComputationMask(maskImage);

    if (meanImage)
        mainFilter->SetMeanImage(meanImage);

    if (varImage)
        mainFilter->SetVarianceImage(varImage);
    
    callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);
    mainFilter->AddObserver(itk::ProgressEvent(), callback);
    
    std::cout << "First run for estimating the scale parameter..." << std::endl;

    mainFilter->Update();

    double scale = mainFilter->GetScale();
    
    std::cout << "\nDone. Estimated scale parameter: " << scale << std::endl;
    
    // Now, transform Rician noise to Gaussian noise using fixed scale parameter
    mainFilter = RiceToGaussianImageFilterType::New();
    mainFilter->SetInput(inputImage);
    mainFilter->SetMaximumNumberOfIterations(args.maxiters);
    mainFilter->SetEpsilon(args.epsilon);
    mainFilter->SetSigma(args.sigma);
    mainFilter->SetScale(scale);
    mainFilter->SetNumberOfWorkUnits(args.nthreads);
    
    if (maskImage)
        mainFilter->SetComputationMask(maskImage);
    
    callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);
    mainFilter->AddObserver(itk::ProgressEvent(), callback);
    
    mainFilter->Update();
    
    std::cout << "\n" << std::endl;

    return mainFilter->GetGaussianImage();
}

template <class ComponentType, unsigned int ImageDimension>
void
transformNoise(itk::ImageIOBase::Pointer imageIO, const arguments &args)
{
    if (ImageDimension == 4)
    {
        typedef itk::Image<ComponentType,4> Image4DType;
        typedef itk::Image<ComponentType,3> Image3DType;
        typename Image4DType::Pointer input = anima::readImage<Image4DType>(args.input);
        std::vector <typename Image3DType::Pointer> inputData = anima::getImagesFromHigherDimensionImage<Image4DType,Image3DType>(input);
        unsigned int nbImages = inputData.size();
        std::vector <typename Image3DType::Pointer> meanData(nbImages), varData(nbImages);
        if (args.mean != "")
            meanData = anima::getImagesFromHigherDimensionImage<Image4DType,Image3DType>(anima::readImage<Image4DType>(args.mean));
        if (args.variance != "")
            varData = anima::getImagesFromHigherDimensionImage<Image4DType,Image3DType>(anima::readImage<Image4DType>(args.variance));
        
        for (unsigned int i = 0;i < nbImages;++i)
            inputData[i] = performTransformation<ComponentType,3>(inputData[i], meanData[i], varData[i], args);
        concatenate<ComponentType>(inputData,args.output);
        return;
    }

    typedef itk::Image<ComponentType,ImageDimension> InputImageType;
    typename InputImageType::Pointer inputImage = anima::readImage<InputImageType>(args.input);
    typename InputImageType::Pointer meanImage, varImage;

    if (args.mean != "")
        meanImage = anima::readImage<InputImageType>(args.mean);

    if (args.variance != "")
        varImage = anima::readImage<InputImageType>(args.variance);

    inputImage = performTransformation<ComponentType,ImageDimension>(inputImage,meanImage,varImage,args);
    
    anima::writeImage<InputImageType>(args.output, inputImage);
}

template <class ComponentType>
void
retrieveNbDimensions(itk::ImageIOBase::Pointer imageIO, const arguments &args)
{
    ANIMA_RETRIEVE_NUMBER_OF_DIMENSIONS(imageIO, ComponentType, transformNoise, imageIO, args)
}

int main(int ac, const char** av)
{
    
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inputArg("i", "input", "Input Rice-corrupted image", true, "", "Input Rice-corrupted image", cmd);
    TCLAP::ValueArg<std::string> outputArg("o", "output", "Output Gaussian-corrupted image", true, "", "Output Gaussian-corrupted image", cmd);
    
    TCLAP::ValueArg<std::string> maskArg("m", "mask", "Segmentation mask", false, "", "Segmentation mask", cmd);
    TCLAP::ValueArg<std::string> meanArg("", "mean", "Mean image", false, "", "Mean image", cmd);
    TCLAP::ValueArg<std::string> varArg("", "variance", "Variance image", false, "", "Variance image", cmd);
    
    TCLAP::ValueArg<double> epsArg("E", "epsilon", "Minimal absolute value difference betweem fixed point iterations (default: 1e-8)", false, 1.0e-8, "Minimal absolute value difference betweem fixed point iterations", cmd);
    TCLAP::ValueArg<double> sigmaArg("S", "sigma", "Gaussian standard deviation for defining neighbor weights (default: 1.0)", false, 1.0, "Gaussian standard deviation for defining neighbor weights", cmd);
    
    TCLAP::ValueArg<unsigned int> maxiterArg("I", "maxiter", "Maximum number of iterations (default: 100)", false, 100, "Maximum number of iterations", cmd);
    TCLAP::ValueArg<unsigned int> nbpArg("T", "nbp", "Number of threads to run on -> default : automatically determine", false, itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(), "Number of threads", cmd);
    
    try
    {
        cmd.parse(ac,av);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    // Retrieve image info
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inputArg.getValue().c_str(),
                                                                           itk::ImageIOFactory::ReadMode);
    if (!imageIO)
    {
        std::cerr << "Itk could not find suitable IO factory for the input" << std::endl;
        return EXIT_FAILURE;
    }
    
    // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
    imageIO->SetFileName(inputArg.getValue());
    imageIO->ReadImageInformation();
    
    arguments args;
    args.input = inputArg.getValue();
    args.output = outputArg.getValue();
    args.mask = maskArg.getValue();
    args.mean = meanArg.getValue();
    args.variance = varArg.getValue();
    args.epsilon = epsArg.getValue();
    args.sigma = sigmaArg.getValue();
    args.maxiters = maxiterArg.getValue();
    args.nthreads = nbpArg.getValue();
    
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
