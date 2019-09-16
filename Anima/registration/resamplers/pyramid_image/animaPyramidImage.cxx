#include <tclap/CmdLine.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <animaLogTensorImageFilter.h>
#include <animaExpTensorImageFilter.h>

#include <animaPyramidImageFilter.h>
#include <animaResampleImageFilter.h>
#include <animaTensorResampleImageFilter.h>

int main(int ac, const char** av)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","input","Input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","Output resampled image",true,"","output image",cmd);

    TCLAP::ValueArg<unsigned int> pyrNumArg("p","pyr-num","number of pyramid levels",false,4,"number of pyramid level",cmd);
    TCLAP::ValueArg<unsigned int> pyrLevelArg("l","pyr-level","pyramid level image required",false,1,"pyramid level",cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("T","numberofthreads","Number of threads to run on (default : all cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(ac,av);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    const    unsigned int    Dimension = 3;
    typedef  float           PixelType;

    typedef itk::Image< PixelType, Dimension >  ImageType;

    typedef itk::ImageFileReader<ImageType> ReaderType;
    typedef itk::ImageFileWriter<ImageType> WriterType;

    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inArg.getValue().c_str(),
                                                                           itk::ImageIOFactory::ReadMode);

    if (!imageIO)
    {
        std::cerr << "Itk could not find suitable IO factory for the input" << std::endl;
        return EXIT_FAILURE;
    }

    // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
    imageIO->SetFileName(inArg.getValue());
    imageIO->ReadImageInformation();

    unsigned int numComponents = imageIO->GetNumberOfComponents();

    bool tensorImage = false;
    if (numComponents == 6)
        tensorImage = true;
    else if (numComponents != 1)
    {
        std::cerr << "Error: not a tensor or scalar image, don't know what to do" << std::endl;
        return(1);
    }

    if (!tensorImage)
    {
        ReaderType::Pointer reader = ReaderType::New();
        reader->SetImageIO(imageIO);
        reader->SetFileName(inArg.getValue());

        reader->Update();

        typedef anima::PyramidImageFilter <ImageType,ImageType> PyramidType;
        typedef anima::ResampleImageFilter<ImageType, ImageType, double> ResampleFilterType;

        PyramidType::Pointer imagePyramid = PyramidType::New();
        imagePyramid->SetInput(reader->GetOutput());
        ResampleFilterType::Pointer resampler = ResampleFilterType::New();
        imagePyramid->SetImageResampler(resampler);
        imagePyramid->SetNumberOfLevels(pyrNumArg.getValue());
        imagePyramid->SetNumberOfWorkUnits(nbpArg.getValue());
        imagePyramid->Update();

        if (pyrLevelArg.getValue() >= imagePyramid->GetNumberOfLevels())
        {
            std::cerr << "Error: unable to compute a pyramid image so far, max is: " << imagePyramid->GetNumberOfLevels() - 1 << std::endl;
            return EXIT_FAILURE;
        }

        ImageType::Pointer outputImage = imagePyramid->GetOutput(pyrLevelArg.getValue());
        outputImage->DisconnectPipeline();

        WriterType::Pointer outWriter = WriterType::New();
        outWriter->SetInput(outputImage);
        outWriter->SetFileName(outArg.getValue());
        outWriter->SetUseCompression(true);
        outWriter->Update();

        return EXIT_SUCCESS;
    }

    typedef itk::VectorImage <PixelType,Dimension> VectorImageType;
    typedef itk::ImageFileReader<VectorImageType> VectorReaderType;
    typedef itk::ImageFileWriter<VectorImageType> VectorWriterType;

    VectorReaderType::Pointer vectorReader = VectorReaderType::New();
    vectorReader->SetImageIO(imageIO);
    vectorReader->Update();

    typedef anima::LogTensorImageFilter <PixelType,Dimension> LogTensorFilterType;
    LogTensorFilterType::Pointer tensorLogger = LogTensorFilterType::New();

    tensorLogger->SetInput(vectorReader->GetOutput());
    tensorLogger->SetScaleNonDiagonal(false);
    tensorLogger->SetNumberOfWorkUnits(nbpArg.getValue());

    tensorLogger->Update();

    typedef anima::PyramidImageFilter <VectorImageType,VectorImageType> VectorPyramidType;
    typedef anima::TensorResampleImageFilter <VectorImageType,double> ResampleFilterType;

    ResampleFilterType::Pointer vectorResampler = ResampleFilterType::New();
    vectorResampler->SetFiniteStrainReorientation(true);

    VectorPyramidType::Pointer vectorImagePyramid = VectorPyramidType::New();
    vectorImagePyramid->SetInput(vectorReader->GetOutput());
    vectorImagePyramid->SetImageResampler(vectorResampler);
    vectorImagePyramid->SetNumberOfLevels(pyrNumArg.getValue());
    vectorImagePyramid->SetNumberOfWorkUnits(nbpArg.getValue());
    vectorImagePyramid->Update();

    if (pyrLevelArg.getValue() >= vectorImagePyramid->GetNumberOfLevels())
    {
        std::cerr << "Error: unable to compute a pyramid image so far, max is: " << vectorImagePyramid->GetNumberOfLevels() - 1 << std::endl;
        return EXIT_FAILURE;
    }

    VectorImageType::Pointer vectorOutputImage = vectorImagePyramid->GetOutput(pyrLevelArg.getValue());
    vectorOutputImage->DisconnectPipeline();

    typedef anima::ExpTensorImageFilter <PixelType,Dimension> ExpTensorFilterType;
    ExpTensorFilterType::Pointer tensorExper = ExpTensorFilterType::New();

    tensorExper->SetInput(vectorOutputImage);
    tensorExper->SetScaleNonDiagonal(false);
    tensorExper->SetNumberOfWorkUnits(nbpArg.getValue());

    tensorExper->Update();

    VectorWriterType::Pointer outVecWriter = VectorWriterType::New();
    outVecWriter->SetInput(tensorExper->GetOutput());
    outVecWriter->SetFileName(outArg.getValue());
    outVecWriter->SetUseCompression(true);
    outVecWriter->Update();

    return EXIT_SUCCESS;
}
