#include <tclap/CmdLine.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <animaTransformSeriesReader.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

#include <animaODFResampleImageFilter.h>

int main(int ac, const char** av)
{
    std::string descriptionMessage;
    descriptionMessage += "Resampler tool to apply a series of transformations to an image of ODF models in Max Descoteaux's basis. ";
    descriptionMessage += "Input transform is an XML file describing all transforms to apply. ";
    descriptionMessage += "Such a file should look like this:\n";
    descriptionMessage += "<TransformationList>\n";
    descriptionMessage += "<Transformation>\n";
    descriptionMessage += "<Type>linear</Type> (it can be svf or dense too)\n";
    descriptionMessage += "<Path>FileName</Path>\n";
    descriptionMessage += "<Inversion>0</Inversion>\n";
    descriptionMessage += "</Transformation>\n";
    descriptionMessage += "...\n";
    descriptionMessage += "</TransformationList>\n\n";
    descriptionMessage += "INRIA / IRISA - VisAGeS/Empenn Team";
    
    TCLAP::CmdLine cmd(descriptionMessage, ' ',ANIMA_VERSION);
	
    TCLAP::ValueArg<std::string> inArg("i","input","Input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> trArg("t","trsf","Transformations XML list",true,"","transformations list",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","Output resampled image",true,"","output image",cmd);
    TCLAP::ValueArg<std::string> geomArg("g","geometry","Geometry image",true,"","geometry image",cmd);
    
    TCLAP::ValueArg<unsigned int> expOrderArg("e","exp-order","Order of field exponentiation approximation (in between 0 and 1, default: 0)",false,0,"exponentiation order",cmd);
    TCLAP::SwitchArg invertArg("I","invert","Invert the transformation series",cmd,false);
    TCLAP::SwitchArg nearestArg("N","nearest","Use nearest neighbor interpolation",cmd,false);
    
    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default: all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
    
    try
    {
        cmd.parse(ac,av);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    const    unsigned int    Dimension = 3;
    typedef  float           PixelType;
    
    typedef itk::VectorImage< PixelType, Dimension >  ImageType;
    typedef anima::TransformSeriesReader <double, Dimension> TransformSeriesReaderType;
    typedef TransformSeriesReaderType::OutputTransformType TransformType;
    
    typedef itk::ImageFileReader <ImageType> ReaderType;
    typedef itk::ImageFileWriter <ImageType> WriterType;

    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(inArg.getValue());
    reader->Update();
    
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(geomArg.getValue().c_str(),
                                                                           itk::ImageIOFactory::ReadMode);

    if (!imageIO)
    {
        std::cerr << "Itk could not find suitable IO factory for the input" << std::endl;
        return EXIT_FAILURE;
    }

    // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
    imageIO->SetFileName(geomArg.getValue());
    imageIO->ReadImageInformation();

    TransformSeriesReaderType *trReader = new TransformSeriesReaderType;
    trReader->SetInput(trArg.getValue());
    trReader->SetInvertTransform(invertArg.isSet());
    trReader->SetExponentiationOrder(expOrderArg.getValue());
    trReader->SetNumberOfThreads(nbpArg.getValue());

    try
    {
        trReader->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }
    
    TransformType::Pointer trsf = trReader->GetOutputTransform();
    
    itk::InterpolateImageFunction <ImageType>::Pointer interpolator;
    
    if (nearestArg.isSet())
        interpolator = itk::NearestNeighborInterpolateImageFunction<ImageType>::New();
    else
        interpolator = anima::VectorModelLinearInterpolateImageFunction<ImageType>::New();
    
    typedef anima::ODFResampleImageFilter <ImageType, double> ResampleFilterType;

    ResampleFilterType::Pointer resample = ResampleFilterType::New();
        
    resample->SetTransform(trsf);
    resample->SetFiniteStrainReorientation(true);
    resample->SetInterpolator(interpolator.GetPointer());
    resample->SetNumberOfThreads(nbpArg.getValue());

    ImageType::DirectionType directionMatrix;
    ImageType::PointType origin;
    ImageType::SpacingType spacing;
    ImageType::RegionType largestRegion;

    for (unsigned int i = 0;i < Dimension;++i)
    {
        origin[i] = imageIO->GetOrigin(i);
        spacing[i] = imageIO->GetSpacing(i);
        largestRegion.SetIndex(i,0);
        largestRegion.SetSize(i,imageIO->GetDimensions(i));

        for (unsigned int j = 0;j < Dimension;++j)
            directionMatrix(i,j) = imageIO->GetDirection(j)[i];
    }

    resample->SetOutputLargestPossibleRegion(largestRegion);
    resample->SetOutputOrigin(origin);
    resample->SetOutputSpacing(spacing);
    resample->SetOutputDirection(directionMatrix);
    
    resample->SetInput(reader->GetOutput());
    
    std::cout << "Applying transform... " << std::flush;
    resample->Update();
    std::cout << "Done..." << std::endl;

    WriterType::Pointer writer = WriterType::New();

    writer->SetUseCompression(true);
    
    writer->SetInput(resample->GetOutput());
    
    writer->SetFileName(outArg.getValue());
    
    writer->Update();
    
    return EXIT_SUCCESS;
}
