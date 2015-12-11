#include <tclap/CmdLine.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <animaTransformSeriesReader.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <animaVectorModelLinearInterpolateImageFunction.h>

#include <animaFiniteStrainTensorResampleImageFilter.h>
#include <animaLogTensorImageFilter.h>
#include <animaExpTensorImageFilter.h>

int main(int ac, const char** av)
{
    std::string descriptionMessage;
    descriptionMessage += "Resampler tool to apply a series of transformations to a tensor image. ";
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
    descriptionMessage += "INRIA / IRISA - VisAGeS Team";

    TCLAP::CmdLine cmd(descriptionMessage, ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","input","Input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> trArg("t","trsf","Transformations XML list",true,"","transformations list",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","Output resampled image",true,"","output image",cmd);
    TCLAP::ValueArg<std::string> geomArg("g","geometry","Geometry image",true,"","geometry image",cmd);

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
        return(1);
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

    ReaderType::Pointer readerGeometry = ReaderType::New();
    readerGeometry->SetFileName(geomArg.getValue());
    readerGeometry->GenerateOutputInformation();

    TransformSeriesReaderType *trReader = new TransformSeriesReaderType;
    trReader->SetInput(trArg.getValue());
    trReader->SetInvertTransform(invertArg.isSet());

    try
    {
        trReader->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return -1;
    }

    TransformType::Pointer trsf = trReader->GetOutputTransform();

    itk::InterpolateImageFunction <ImageType>::Pointer interpolator;

    if (nearestArg.isSet())
        interpolator = itk::NearestNeighborInterpolateImageFunction<ImageType>::New();
    else
        interpolator = anima::VectorModelLinearInterpolateImageFunction<ImageType>::New();

    typedef anima::FiniteStrainTensorResampleImageFilter <ImageType, double> ResampleFilterType;

    ResampleFilterType::Pointer resample = ResampleFilterType::New();

    resample->SetTransform(trsf);
    resample->SetInterpolator(interpolator.GetPointer());
    resample->SetNumberOfThreads(nbpArg.getValue());

    resample->SetOutputLargestPossibleRegion(readerGeometry->GetOutput()->GetLargestPossibleRegion());
    resample->SetOutputOrigin(readerGeometry->GetOutput()->GetOrigin());
    resample->SetOutputSpacing(readerGeometry->GetOutput()->GetSpacing());
    resample->SetOutputDirection(readerGeometry->GetOutput()->GetDirection());

    typedef anima::LogTensorImageFilter <PixelType,Dimension> LogTensorFilterType;
    LogTensorFilterType::Pointer tensorLogger = LogTensorFilterType::New();

    tensorLogger->SetInput(reader->GetOutput());
    tensorLogger->SetScaleNonDiagonal(true);
    tensorLogger->SetNumberOfThreads(nbpArg.getValue());

    std::cout << "Logging input... " << std::flush;
    tensorLogger->Update();
    std::cout << "Done..." << std::endl;

    ImageType::Pointer tmpImage = tensorLogger->GetOutput();
    tmpImage->DisconnectPipeline();

    resample->SetInput(tmpImage);

    std::cout << "Applying transform... " << std::flush;
    resample->Update();
    std::cout << "Done..." << std::endl;

    tmpImage = resample->GetOutput();
    tmpImage->DisconnectPipeline();

    typedef anima::ExpTensorImageFilter <PixelType,Dimension> ExpTensorFilterType;
    ExpTensorFilterType::Pointer tensorExper = ExpTensorFilterType::New();

    tensorExper->SetInput(tmpImage);
    tensorExper->SetScaleNonDiagonal(true);
    tensorExper->SetNumberOfThreads(nbpArg.getValue());

    std::cout << "Exping output... " << std::flush;
    tensorExper->Update();
    std::cout << "Done..." << std::endl;

    WriterType::Pointer writer = WriterType::New();

    writer->SetUseCompression(true);

    tmpImage = tensorExper->GetOutput();
    tmpImage->DisconnectPipeline();

    writer->SetInput(tmpImage);

    writer->SetFileName(outArg.getValue());

    writer->Update();

    return 0;
}
