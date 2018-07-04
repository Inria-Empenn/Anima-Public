#include <animaRiceToGaussianImageFilter.h>
#include <animaReadWriteFunctions.h>
#include <animaGradientFileReader.h>

#include <tclap/CmdLine.h>

#include <itkSubtractImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkNaryAddImageFilter.h>

//Update progression of the process
void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
    std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
}

std::string to_string(const unsigned int n, const unsigned int width)
{
    std::ostringstream stm;
    
    if (width == 0)
        stm << n;
    else
        stm << std::setw(width) << std::setfill('0') << n;
    
    return stm.str();
}

int main(int ac, const char** av)
{
    
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inputArg("i", "input", "Input Rice-corrupted image", true, "", "Input Rice-corrupted image", cmd);
    TCLAP::ValueArg<std::string> outputArg("o", "output", "Output Gaussian-corrupted image", true, "", "Output Gaussian-corrupted image", cmd);
    TCLAP::ValueArg<std::string> bvecArg("g","grad","Input gradients",true,"","Input gradients",cmd);
    TCLAP::ValueArg<std::string> bvalArg("b","bval","Input b-values",true,"","Input b-values",cmd);
    
    TCLAP::ValueArg<std::string> scaleArg("s", "scale", "Output scale image", false, "", "Output scale image", cmd);
    TCLAP::ValueArg<std::string> maskArg("m", "mask", "Segmentation mask", false, "", "Segmentation mask", cmd);
    
    TCLAP::ValueArg<double> epsArg("E", "epsilon", "Minimal absolute value difference betweem fixed point iterations (default: 1e-8)", false, 1.0e-8, "Minimal absolute value difference betweem fixed point iterations", cmd);
    TCLAP::ValueArg<double> sigmaArg("S", "sigma", "Gaussian standard deviation for defining neighbor weights (default: 1.0)", false, 1.0, "Gaussian standard deviation for defining neighbor weights", cmd);
    
    TCLAP::ValueArg<unsigned int> maxiterArg("I", "maxiter", "Maximum number of iterations (default: 100)", false, 100, "Maximum number of iterations", cmd);
    TCLAP::ValueArg<unsigned int> nbpArg("T", "nbp", "Number of threads to run on -> default : automatically determine", false, itk::MultiThreader::GetGlobalDefaultNumberOfThreads(), "Number of threads", cmd);
    
    try
    {
        cmd.parse(ac,av);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }
    
    const unsigned int ImageDimension = 3;
    
    typedef anima::GradientFileReader <vnl_vector_fixed<double,ImageDimension>,double> GFReaderType;
    GFReaderType gfReader;
    gfReader.SetGradientFileName(bvecArg.getValue());
    gfReader.SetBValueBaseString(bvalArg.getValue());
    gfReader.Update();
    
    GFReaderType::GradientVectorType directions = gfReader.GetGradients();
    GFReaderType::BValueVectorType bvalues = gfReader.GetBValues();
    
    typedef itk::Image<float,4> Image4DType;
    typedef anima::RiceToGaussianImageFilter<ImageDimension> RiceToGaussianImageFilterType;
    typedef RiceToGaussianImageFilterType::InputImageType InputImageType;
    typedef RiceToGaussianImageFilterType::OutputImageType OutputImageType;
    typedef RiceToGaussianImageFilterType::MaskImageType MaskImageType;
    typedef itk::SubtractImageFilter<InputImageType,InputImageType,OutputImageType> SubtractFilterType;
    typedef itk::MultiplyImageFilter<InputImageType,InputImageType,OutputImageType> MultiplyImageFilterType;
    typedef itk::NaryAddImageFilter<InputImageType,OutputImageType> NaryAddImageFilterType;
    
    Image4DType::Pointer input = anima::readImage<Image4DType>(inputArg.getValue());
    std::vector <InputImageType::Pointer> inputData;
    inputData = anima::getImagesFromHigherDimensionImage<Image4DType,InputImageType>(input);
    
    unsigned int numberOfImages = inputData.size();
    unsigned int numberOfB0Images = 0;
    unsigned int indexOfFirstB0Image = 0;
    
    for (unsigned int i = 0;i < numberOfImages;++i)
    {
        if (bvalues[i] <= 10)
        {
            if (numberOfB0Images == 0)
                indexOfFirstB0Image = i;
            ++numberOfB0Images;
        }
    }
    
    std::cout << "Number of B0 images: " << numberOfB0Images << std::endl;
    
    itk::CStyleCommand::Pointer callback;
    RiceToGaussianImageFilterType::Pointer mainFilter;
    
    // Estimate scale parameter from B0 image(s)
    if (numberOfB0Images == 1)
    {
        mainFilter = RiceToGaussianImageFilterType::New();
        mainFilter->SetInput(inputData[indexOfFirstB0Image]);
        mainFilter->SetMaximumNumberOfIterations(maxiterArg.getValue());
        mainFilter->SetEpsilon(epsArg.getValue());
        mainFilter->SetSigma(sigmaArg.getValue());
        mainFilter->SetNumberOfThreads(nbpArg.getValue());
        
        if (maskArg.getValue() != "")
            mainFilter->SetComputationMask(anima::readImage<MaskImageType>(maskArg.getValue()));
        
        callback = itk::CStyleCommand::New();
        callback->SetCallback(eventCallback);
        mainFilter->AddObserver(itk::ProgressEvent(), callback);
        
        std::cout << "First run of Gaussianizer on the B0 image for estimating the scale parameter..." << std::endl;
    }
    else
    {
        NaryAddImageFilterType::Pointer addFilter = NaryAddImageFilterType::New();
        unsigned int pos = 0;
        for (unsigned int i = 0;i < numberOfImages;++i)
        {
            if (bvalues[i] > 10)
                continue;
            
            addFilter->SetInput(pos, inputData[i]);
            ++pos;
        }
        
        addFilter->SetNumberOfThreads(nbpArg.getValue());
        addFilter->Update();
        
        MultiplyImageFilterType::Pointer mulFilter = MultiplyImageFilterType::New();
        mulFilter->SetInput1(addFilter->GetOutput());
        mulFilter->SetConstant2(1.0 / (double)numberOfB0Images);
        mulFilter->SetNumberOfThreads(nbpArg.getValue());
        mulFilter->Update();
        
        InputImageType::Pointer meanImage = mulFilter->GetOutput();
        
        addFilter = NaryAddImageFilterType::New();
        pos = 0;
        for (unsigned int i = 0;i < numberOfImages;++i)
        {
            if (bvalues[i] > 10)
                continue;
            
            SubtractFilterType::Pointer subFilter = SubtractFilterType::New();
            subFilter->SetInput1(inputData[i]);
            subFilter->SetInput2(meanImage);
            subFilter->SetNumberOfThreads(nbpArg.getValue());
            subFilter->Update();
            
            mulFilter = MultiplyImageFilterType::New();
            mulFilter->SetInput1(subFilter->GetOutput());
            mulFilter->SetInput2(subFilter->GetOutput());
            mulFilter->SetNumberOfThreads(nbpArg.getValue());
            mulFilter->Update();
            
            addFilter->SetInput(pos, mulFilter->GetOutput());
            ++pos;
        }
        
        addFilter->SetNumberOfThreads(nbpArg.getValue());
        addFilter->Update();
        
        mulFilter = MultiplyImageFilterType::New();
        mulFilter->SetInput1(addFilter->GetOutput());
        mulFilter->SetConstant2(1.0 / (numberOfB0Images - 1.0));
        mulFilter->SetNumberOfThreads(nbpArg.getValue());
        mulFilter->Update();
        
        InputImageType::Pointer varianceImage = mulFilter->GetOutput();
        
        mainFilter = RiceToGaussianImageFilterType::New();
        mainFilter->SetInput(inputData[indexOfFirstB0Image]);
        mainFilter->SetMaximumNumberOfIterations(maxiterArg.getValue());
        mainFilter->SetEpsilon(epsArg.getValue());
        mainFilter->SetSigma(sigmaArg.getValue());
        mainFilter->SetMeanImage(meanImage);
        mainFilter->SetVarianceImage(varianceImage);
        mainFilter->SetNumberOfThreads(nbpArg.getValue());
        
        if (maskArg.getValue() != "")
            mainFilter->SetComputationMask(anima::readImage<MaskImageType>(maskArg.getValue()));
        
        callback = itk::CStyleCommand::New();
        callback->SetCallback(eventCallback);
        mainFilter->AddObserver(itk::ProgressEvent(), callback);
        
        std::cout << "First run of Gaussianizer using the mean and variance of the B0 images for estimating the scale parameter..." << std::endl;
    }
    
    try
    {
        mainFilter->Update();
    }
    catch (itk::ExceptionObject & err)
    {
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }
    
    double scale = mainFilter->GetScale();
    
    if (scaleArg.getValue() != "")
        anima::writeImage<OutputImageType>(scaleArg.getValue(), mainFilter->GetScaleImage());
    
    std::cout << "\nDone. Estimated scale parameter: " << scale << std::endl;
    
    // Transform Rice to Gaussian data
    
    for (unsigned int i = 0;i < numberOfImages;++i)
    {
        std::cout << "\nRun Gaussianizer for Image "<< i+1 << " with fixed scale parameter..." << std::endl;
        
        mainFilter = RiceToGaussianImageFilterType::New();
        mainFilter->SetInput(inputData[i]);
        mainFilter->SetMaximumNumberOfIterations(maxiterArg.getValue());
        mainFilter->SetEpsilon(epsArg.getValue());
        mainFilter->SetSigma(sigmaArg.getValue());
        mainFilter->SetScale(scale);
        mainFilter->SetNumberOfThreads(nbpArg.getValue());
        
        if (maskArg.getValue() != "")
            mainFilter->SetComputationMask(anima::readImage<MaskImageType>(maskArg.getValue()));
        
        callback = itk::CStyleCommand::New();
        callback->SetCallback(eventCallback);
        mainFilter->AddObserver(itk::ProgressEvent(), callback);
        
        try
        {
            mainFilter->Update();
        }
        catch (itk::ExceptionObject & err)
        {
            std::cerr << err << std::endl;
            return EXIT_FAILURE;
        }
        
        anima::writeImage<OutputImageType>(outputArg.getValue() + to_string(i + 1, 3) + ".nrrd", mainFilter->GetGaussianImage());
    }
    
    return EXIT_SUCCESS;
}
