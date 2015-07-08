#include <tclap/CmdLine.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkImageRegionIterator.h>

#include <itkAddImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkMultiplyImageFilter.h>

template <class ImageType>
void computeArithmetic (std::string &inStr, std::string &outStr, std::string &multImStr, std::string &divImStr, std::string &addImStr,
                        std::string &subImStr, double multConstant, double divideConstant, double addConstant,
                        unsigned int nThreads)
{
    typedef itk::ImageFileReader <ImageType> ReaderType;
    typedef itk::ImageFileWriter <ImageType> WriterType;
    typedef itk::AddImageFilter <ImageType,ImageType,ImageType> AddFilterType;
    typedef itk::AddImageFilter <ImageType,itk::Image <double, ImageType::ImageDimension>,ImageType> AddConstantFilterType;
    typedef itk::SubtractImageFilter <ImageType,ImageType,ImageType> SubtractFilterType;
    
    typedef itk::MultiplyImageFilter <ImageType,itk::Image <double, ImageType::ImageDimension>,ImageType> MultiplyFilterType;
    typedef itk::DivideImageFilter <ImageType,itk::Image <double, ImageType::ImageDimension>,ImageType> DivideFilterType;
    
    typename ReaderType::Pointer inputReader = ReaderType::New();
    inputReader->SetFileName(inStr);
    inputReader->Update();
    
    typename ImageType::Pointer currentImage = inputReader->GetOutput();
    currentImage->DisconnectPipeline();
    
    if (multImStr != "")
    {
        typedef itk::ImageFileReader < itk::Image <double, ImageType::ImageDimension> > MultiplyReaderType;
        typename MultiplyReaderType::Pointer multImReader = MultiplyReaderType::New();
        multImReader->SetFileName(multImStr);
        multImReader->Update();
        
        typename MultiplyFilterType::Pointer multFilter = MultiplyFilterType::New();
        multFilter->SetInput1(currentImage);
        multFilter->SetInput2(multImReader->GetOutput());
        multFilter->SetNumberOfThreads(nThreads);
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
        multFilter->SetNumberOfThreads(nThreads);
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
        divFilter->SetNumberOfThreads(nThreads);
        divFilter->InPlaceOn();
        
        divFilter->Update();
        currentImage = divFilter->GetOutput();
        currentImage->DisconnectPipeline();
    }
    
    if (divImStr != "")
    {
        typedef itk::ImageFileReader < itk::Image <double, ImageType::ImageDimension> > DivideReaderType;
        typename DivideReaderType::Pointer divImReader = DivideReaderType::New();
        divImReader->SetFileName(divImStr);
        divImReader->Update();

        typename DivideFilterType::Pointer divFilter = DivideFilterType::New();
        divFilter->SetInput1(currentImage);
        divFilter->SetInput2(divImReader->GetOutput());
        divFilter->SetNumberOfThreads(nThreads);
        divFilter->InPlaceOn();
        
        divFilter->Update();
        currentImage = divFilter->GetOutput();
        currentImage->DisconnectPipeline();
        
        itk::ImageRegionIterator < itk::Image <double, ImageType::ImageDimension> > divImItr (divImReader->GetOutput(),divImReader->GetOutput()->GetLargestPossibleRegion());
        itk::ImageRegionIterator <ImageType> currentImageItr (currentImage,currentImage->GetLargestPossibleRegion());
        
        while (!currentImageItr.IsAtEnd())
        {
            if (divImItr.Get() == 0)
                currentImageItr.Set(divImItr.Get() * currentImageItr.Get());
            
            ++currentImageItr;
            ++divImItr;
        }
    }
    
    if (addImStr != "")
    {
        typename ReaderType::Pointer addImReader = ReaderType::New();
        addImReader->SetFileName(addImStr);
        addImReader->Update();
        
        typename AddFilterType::Pointer addFilter = AddFilterType::New();
        addFilter->SetInput1(currentImage);
        addFilter->SetInput2(addImReader->GetOutput());
        addFilter->SetNumberOfThreads(nThreads);
        addFilter->InPlaceOn();
        
        addFilter->Update();
        currentImage = addFilter->GetOutput();
        currentImage->DisconnectPipeline();
    }
    
    if (subImStr != "")
    {
        typename ReaderType::Pointer subImReader = ReaderType::New();
        subImReader->SetFileName(subImStr);
        subImReader->Update();
        
        typename SubtractFilterType::Pointer subFilter = SubtractFilterType::New();
        subFilter->SetInput1(currentImage);
        subFilter->SetInput2(subImReader->GetOutput());
        subFilter->SetNumberOfThreads(nThreads);
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
        addFilter->SetNumberOfThreads(nThreads);
        addFilter->InPlaceOn();
        
        addFilter->Update();
        currentImage = addFilter->GetOutput();
        currentImage->DisconnectPipeline();
    }
    
    // Finally write the result
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput(currentImage);
    writer->SetFileName(outStr);
    writer->SetUseCompression(true);
    writer->Update();
}

int main(int argc, char **argv)
{
    std::string descriptionMessage = "Performs very basic mathematical operations on images: performs (I * m * M) / (D * d) + A + a - s\n";
    descriptionMessage += "This software has known limitations: you might have to use it several times in a row to perform the operation you want,\n";
    descriptionMessage += "it requires the divide and multiply images to be scalar, and the add and subtract images to be of the same format as the input (although this is not verified).\n";
    descriptionMessage += "INRIA / IRISA - VisAGeS Team";
    
    TCLAP::CmdLine cmd(descriptionMessage, ' ',"1.0");
    
    TCLAP::ValueArg<std::string> inArg("i","input","Input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","output image",true,"","output image",cmd);
    
    TCLAP::ValueArg<std::string> multImArg("m","mult-image","multiply image",false,"","multiply image",cmd);
    TCLAP::ValueArg<std::string> divImArg("d","div-image","divide image",false,"","divide image",cmd);
    TCLAP::ValueArg<std::string> addImArg("a","add-image","add image",false,"","add image",cmd);
    TCLAP::ValueArg<std::string> subtractImArg("s","sub-image","subtract image",false,"","subtract image",cmd);

    TCLAP::ValueArg<double> multiplyConstantArg("M","multiply-constant","multiply constant value",false,1.0,"multiply constant value",cmd);
    TCLAP::ValueArg<double> divideConstantArg("D","divide-constant","divide constant value",false,1.0,"divide constant value",cmd);
    TCLAP::ValueArg<double> addConstantArg("A","add-constant","add constant value",false,0.0,"add constant value",cmd);
	
    TCLAP::ValueArg<unsigned int> nbpArg("T","numberofthreads","Number of threads to run on (default : all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
    
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
    typedef itk::ImageFileReader <ImageType> ImageReaderType;
	typedef itk::VectorImage <float,3> VectorImageType;

    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(inArg.getValue());
    reader->GenerateOutputInformation();
    
    bool vectorImage = (reader->GetImageIO()->GetNumberOfComponents() > 1);
    bool fourDimensionalImage = (reader->GetImageIO()->GetNumberOfDimensions() == 4);

    if (vectorImage)
        computeArithmetic<VectorImageType>(inArg.getValue(),outArg.getValue(),multImArg.getValue(),divImArg.getValue(),addImArg.getValue(),
                                           subtractImArg.getValue(),multiplyConstantArg.getValue(),divideConstantArg.getValue(),
                                           addConstantArg.getValue(),nbpArg.getValue());
    else if (fourDimensionalImage)
        computeArithmetic<Image4DType>(inArg.getValue(),outArg.getValue(),multImArg.getValue(),divImArg.getValue(),addImArg.getValue(),
                                       subtractImArg.getValue(),multiplyConstantArg.getValue(),divideConstantArg.getValue(),
                                       addConstantArg.getValue(),nbpArg.getValue());
    else
        computeArithmetic<ImageType>(inArg.getValue(),outArg.getValue(),multImArg.getValue(),divImArg.getValue(),addImArg.getValue(),
                                     subtractImArg.getValue(),multiplyConstantArg.getValue(),divideConstantArg.getValue(),
                                     addConstantArg.getValue(),nbpArg.getValue());
    
    return 0;
}
