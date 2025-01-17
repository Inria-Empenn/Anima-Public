#pragma once

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkExtractImageFilter.h>

namespace anima
{

template <class ImageType>
typename itk::SmartPointer<ImageType>
readImage(std::string filename)
{
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(filename);

    typename itk::SmartPointer<ImageType> img;

    reader->Update();

    img = reader->GetOutput();
    return img;
}

template <class OutputImageType>
void
writeImage(std::string filename, OutputImageType* img)
{
    typedef itk::ImageFileWriter<OutputImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetUseCompression(true);
    writer->SetFileName(filename);
    writer->SetInput(img);

    writer->Update();
}

//! Get a vector of input images from a higher dimensional image
template <class InputImageType, class OutputImageType>
std::vector < itk::SmartPointer <OutputImageType> >
getImagesFromHigherDimensionImage(InputImageType *inputImage)
{
    unsigned int highDimImage = InputImageType::ImageDimension;
    unsigned int lowerDimImage = OutputImageType::ImageDimension;

    if (highDimImage != lowerDimImage + 1)
        throw itk::ExceptionObject(__FILE__, __LINE__, "Trying to divide an image that doesn't have one more dimension",ITK_LOCATION);

    unsigned int ndim = inputImage->GetLargestPossibleRegion().GetSize()[lowerDimImage];

    typename InputImageType::RegionType largeRegion = inputImage->GetLargestPossibleRegion();
    typename InputImageType::RegionType smallRegion = largeRegion;
    typedef itk::ExtractImageFilter <InputImageType, OutputImageType> ExtractFilterType;

    std::vector < itk::SmartPointer <OutputImageType> > outputData;
    smallRegion.SetSize(lowerDimImage,0);
    for (unsigned int i = 0;i < ndim;++i)
    {
        smallRegion.SetIndex(lowerDimImage,i + largeRegion.GetIndex(lowerDimImage));

        typename ExtractFilterType::Pointer extractor = ExtractFilterType::New();
        extractor->SetInput(inputImage);
        extractor->SetExtractionRegion(smallRegion);
        extractor->SetDirectionCollapseToGuess();

        extractor->Update();

        outputData.push_back(extractor->GetOutput());
        outputData[i]->DisconnectPipeline();
    }

    return outputData;
}

//! Set inputs of an image to image filter from a file name containing either a list of files or a higher dimensional image
template <class InputImageType, class ImageFilterType>
unsigned int
setMultipleImageFilterInputsFromFileName(std::string &fileName,
                                         ImageFilterType *filter)
{
    typedef itk::Image <typename InputImageType::PixelType, InputImageType::ImageDimension + 1> HigherDimImageType;
    typedef itk::ImageFileReader <HigherDimImageType> HigherDimImageReaderType;
    typedef itk::ImageFileReader < InputImageType > ImageReaderType;

    unsigned int nbPats = 0;

    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(fileName.c_str(), itk::IOFileModeEnum::ReadMode);

    if( !imageIO ) // file list
    {
        std::ifstream fileIn(fileName.c_str());
        if (!fileIn.is_open())
        {
            std::string errStr = "Unable to read file: ";
            errStr += fileName;

            throw itk::ExceptionObject(__FILE__, __LINE__,errStr,ITK_LOCATION);
        }

        typename ImageReaderType::Pointer imageReader;

        while (!fileIn.eof())
        {
            char tmpStr[2048];
            fileIn.getline(tmpStr,2048);

            if (strcmp(tmpStr,"") == 0)
                continue;

            std::cout << "Loading image " << nbPats << " " << tmpStr << "..." << std::endl;
            imageReader = ImageReaderType::New();
            imageReader->SetFileName(tmpStr);
            imageReader->Update();

            filter->SetInput(nbPats,imageReader->GetOutput());

            nbPats++;
        }

        fileIn.close();
    }
    else // N+1D image tentative
    {
        // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
        imageIO->SetFileName(fileName);
        imageIO->ReadImageInformation();

        // Do we have a (N+1)D image ? If not, exiting, otherwise considering the last dimension as each input
        unsigned int ndim = imageIO->GetNumberOfDimensions();

        if (ndim == InputImageType::ImageDimension)
        {
            filter->SetInput(anima::readImage <InputImageType> (fileName));
            return 1;
        }

        if (ndim != InputImageType::ImageDimension + 1)
        {
            std::string errStr = "Unable to read file: ";
            errStr += fileName;

            throw itk::ExceptionObject(__FILE__, __LINE__,errStr,ITK_LOCATION);
        }

        typename HigherDimImageReaderType::Pointer imageReader = HigherDimImageReaderType::New();
        imageReader->SetImageIO(imageIO);
        imageReader->SetFileName(fileName);
        imageReader->Update();

        std::vector <typename InputImageType::Pointer> inputData;
        inputData = anima::getImagesFromHigherDimensionImage<HigherDimImageType,InputImageType>(imageReader->GetOutput());

        for (unsigned int i = 0;i < inputData.size();++i)
            filter->SetInput(i,inputData[i]);

        nbPats = inputData.size();
    }

    return nbPats;
}

}// end of namespace anima
