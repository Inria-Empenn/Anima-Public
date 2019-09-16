#include <tclap/CmdLine.h>

#include <animaReadWriteFunctions.h>
#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkImageRegionIterator.h>

#include <itkAddImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkPowImageFilter.h>

template <class ImageType>
void computeArithmetic (std::string &inStr, std::string &outStr, std::string &multImStr, std::string &divImStr, std::string &addImStr,
                        std::string &subImStr, double multConstant, double divideConstant, double addConstant, double powConstant,
                        unsigned int nThreads)
{
    typedef itk::Image <double, ImageType::ImageDimension> WorkImageType;
    typedef itk::AddImageFilter <ImageType,ImageType,ImageType> AddFilterType;
    typedef itk::AddImageFilter <ImageType,WorkImageType,ImageType> AddConstantFilterType;
    typedef itk::SubtractImageFilter <ImageType,ImageType,ImageType> SubtractFilterType;
    
    typedef itk::MultiplyImageFilter <ImageType,WorkImageType,ImageType> MultiplyFilterType;
    typedef itk::DivideImageFilter <ImageType,WorkImageType,ImageType> DivideFilterType;
        
    typename ImageType::Pointer currentImage = anima::readImage <ImageType> (inStr);
    currentImage->DisconnectPipeline();
    
    if (multImStr != "")
    {
        typedef itk::ImageFileReader <WorkImageType> MultiplyReaderType;
        typename MultiplyReaderType::Pointer multImReader = MultiplyReaderType::New();
        multImReader->SetFileName(multImStr);
        multImReader->Update();
        
        typename MultiplyFilterType::Pointer multFilter = MultiplyFilterType::New();
        multFilter->SetInput1(currentImage);
        multFilter->SetInput2(multImReader->GetOutput());
        multFilter->SetNumberOfWorkUnits(nThreads);
        multFilter->InPlaceOn();
        
        multFilter->Update();
        currentImage = multFilter->GetOutput();
        currentImage->DisconnectPipeline();
    }
    
    if (multConstant != 1.0)
    {
        typename MultiplyFilterType::Pointer multFilter = MultiplyFilterType::New();
        multFilter->SetInput1(currentImage);
        multFilter->SetConstant(multConstant);
        multFilter->SetNumberOfWorkUnits(nThreads);
        multFilter->InPlaceOn();
        
        multFilter->Update();
        currentImage = multFilter->GetOutput();
        currentImage->DisconnectPipeline();
    }
    
    if ((divideConstant != 0.0)&&(divideConstant != 1.0))
    {
        typename DivideFilterType::Pointer divFilter = DivideFilterType::New();
        divFilter->SetInput1(currentImage);
        divFilter->SetConstant(divideConstant);
        divFilter->SetNumberOfWorkUnits(nThreads);
        divFilter->InPlaceOn();
        
        divFilter->Update();
        currentImage = divFilter->GetOutput();
        currentImage->DisconnectPipeline();
    }
    
    if (divImStr != "")
    {
        typedef itk::ImageFileReader <WorkImageType> DivideReaderType;
        typename DivideReaderType::Pointer divImReader = DivideReaderType::New();
        divImReader->SetFileName(divImStr);
        divImReader->Update();

        typename DivideFilterType::Pointer divFilter = DivideFilterType::New();
        divFilter->SetInput1(currentImage);
        divFilter->SetInput2(divImReader->GetOutput());
        divFilter->SetNumberOfWorkUnits(nThreads);
        divFilter->InPlaceOn();
        
        divFilter->Update();
        currentImage = divFilter->GetOutput();
        currentImage->DisconnectPipeline();
        
        itk::ImageRegionIterator <WorkImageType> divImItr (divImReader->GetOutput(),divImReader->GetOutput()->GetLargestPossibleRegion());
        itk::ImageRegionIterator <ImageType> currentImageItr (currentImage,currentImage->GetLargestPossibleRegion());
        
        while (!currentImageItr.IsAtEnd())
        {
            // Multiplication by 0.0 since current image may be a vector
            if (itk::Math::AlmostEquals(divImItr.Get(), itk::NumericTraits<typename WorkImageType::PixelType>::ZeroValue()))
                currentImageItr.Set(currentImageItr.Get() * 0.0);
            
            ++currentImageItr;
            ++divImItr;
        }
    }
    
    if (addImStr != "")
    {
        typename ImageType::Pointer addedImage = anima::readImage <ImageType> (addImStr);
        
        typename AddFilterType::Pointer addFilter = AddFilterType::New();
        addFilter->SetInput1(currentImage);
        addFilter->SetInput2(addedImage);
        addFilter->SetNumberOfWorkUnits(nThreads);
        addFilter->InPlaceOn();
        
        addFilter->Update();
        currentImage = addFilter->GetOutput();
        currentImage->DisconnectPipeline();
    }
    
    if (subImStr != "")
    {
        typename ImageType::Pointer subtractedImage = anima::readImage <ImageType> (subImStr);
        
        typename SubtractFilterType::Pointer subFilter = SubtractFilterType::New();
        subFilter->SetInput1(currentImage);
        subFilter->SetInput2(subtractedImage);
        subFilter->SetNumberOfWorkUnits(nThreads);
        subFilter->InPlaceOn();
        
        subFilter->Update();
        currentImage = subFilter->GetOutput();
        currentImage->DisconnectPipeline();
    }
    
    if (addConstant != 0.0)
    {
        typename AddConstantFilterType::Pointer addFilter = AddConstantFilterType::New();
        addFilter->SetInput1(currentImage);
        addFilter->SetConstant(addConstant);
        addFilter->SetNumberOfWorkUnits(nThreads);
        addFilter->InPlaceOn();
        
        addFilter->Update();
        currentImage = addFilter->GetOutput();
        currentImage->DisconnectPipeline();
    }
    
    if (powConstant != 1.0)
    {
        typedef itk::Image <typename ImageType::IOPixelType, ImageType::ImageDimension> ScalarImageType;
        typedef itk::PowImageFilter <ScalarImageType, WorkImageType, ScalarImageType> PowConstantFilterType;

        ScalarImageType *castImage = dynamic_cast <ScalarImageType *> (currentImage.GetPointer());
        if (castImage)
        {
            typename PowConstantFilterType::Pointer powFilter = PowConstantFilterType::New();
            powFilter->SetInput1(castImage);
            powFilter->SetConstant(powConstant);
            powFilter->SetNumberOfWorkUnits(nThreads);
            powFilter->InPlaceOn();

            powFilter->Update();

            currentImage = dynamic_cast <ImageType *> (powFilter->GetOutput());
            currentImage->DisconnectPipeline();
        }
    }

    // Finally write the result
    anima::writeImage <ImageType> (outStr,currentImage);
}

int main(int argc, char **argv)
{
    std::string descriptionMessage = "Performs very basic mathematical operations on images: performs ( (I * m * M) / (D * d) + A + a - s )^P \n";
    descriptionMessage += "This software has known limitations: you might have to use it several times in a row to perform the operation you want,\n";
    descriptionMessage += "it requires the divide and multiply images to be scalar, and the add and subtract images to be of the same format as the input (although this is not verified).\n";
    descriptionMessage += "INRIA / IRISA - VisAGeS/Empenn Team";
    
    TCLAP::CmdLine cmd(descriptionMessage, ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inArg("i","input","Input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","output image",true,"","output image",cmd);
    
    TCLAP::ValueArg<std::string> multImArg("m","mult-image","multiply image",false,"","multiply image",cmd);
    TCLAP::ValueArg<std::string> divImArg("d","div-image","divide image",false,"","divide image",cmd);
    TCLAP::ValueArg<std::string> addImArg("a","add-image","add image",false,"","add image",cmd);
    TCLAP::ValueArg<std::string> subtractImArg("s","sub-image","subtract image",false,"","subtract image",cmd);

    TCLAP::ValueArg<double> multiplyConstantArg("M","multiply-constant","multiply constant value",false,1.0,"multiply constant value",cmd);
    TCLAP::ValueArg<double> divideConstantArg("D","divide-constant","divide constant value",false,1.0,"divide constant value",cmd);
    TCLAP::ValueArg<double> addConstantArg("A","add-constant","add constant value",false,0.0,"add constant value",cmd);
    TCLAP::ValueArg<double> powArg("P","pow-constant","power constant value (only for scalar images)",false,1.0,"power constant value",cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("T","numberofthreads","Number of threads to run on (default : all cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
    
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }
	
    typedef itk::Image <float,3> ImageType;
    typedef itk::Image <float,4> Image4DType;
	typedef itk::VectorImage <float,3> VectorImageType;

    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inArg.getValue().c_str(),
                                                                           itk::ImageIOFactory::ReadMode);

    if (!imageIO)
    {
        std::cerr << "Unable to read input image " << inArg.getValue() << std::endl;
        return EXIT_FAILURE;
    }

    imageIO->SetFileName(inArg.getValue());
    imageIO->ReadImageInformation();
    
    bool vectorImage = (imageIO->GetNumberOfComponents() > 1);
    bool fourDimensionalImage = (imageIO->GetNumberOfDimensions() == 4);

    if (vectorImage)
        computeArithmetic<VectorImageType>(inArg.getValue(),outArg.getValue(),multImArg.getValue(),divImArg.getValue(),addImArg.getValue(),
                                           subtractImArg.getValue(),multiplyConstantArg.getValue(),divideConstantArg.getValue(),
                                           addConstantArg.getValue(),powArg.getValue(),nbpArg.getValue());
    else if (fourDimensionalImage)
        computeArithmetic<Image4DType>(inArg.getValue(),outArg.getValue(),multImArg.getValue(),divImArg.getValue(),addImArg.getValue(),
                                       subtractImArg.getValue(),multiplyConstantArg.getValue(),divideConstantArg.getValue(),
                                       addConstantArg.getValue(),powArg.getValue(),nbpArg.getValue());
    else
        computeArithmetic<ImageType>(inArg.getValue(),outArg.getValue(),multImArg.getValue(),divImArg.getValue(),addImArg.getValue(),
                                     subtractImArg.getValue(),multiplyConstantArg.getValue(),divideConstantArg.getValue(),
                                     addConstantArg.getValue(),powArg.getValue(),nbpArg.getValue());
    
    return 0;
}
