#include <tclap/CmdLine.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>

#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkWindowedSincInterpolateImageFunction.h>
#include <itkConstantBoundaryCondition.h>

#include <animaFasterLinearInterpolateImageFunction.h>
#include <animaResampleImageFilter.h>
#include <animaTransformSeriesReader.h>
#include <animaReadWriteFunctions.h>
#include <animaRetrieveImageTypeMacros.h>


struct arguments
{
    bool invert;
    unsigned int pthread;
    std::string input, output, geometry, transfo, interpolation;

};

template <class ImageType>
void
applyVectorTransfo(itk::ImageIOBase::Pointer geometryImageIO, const arguments &args)
{
    typedef itk::VectorImage<float, ImageType::ImageDimension> OuputType;

    typedef anima::TransformSeriesReader <double, ImageType::ImageDimension> TransformSeriesReaderType;
    typedef typename TransformSeriesReaderType::OutputTransformType TransformType;
    TransformSeriesReaderType *trReader = new TransformSeriesReaderType;
    trReader->SetInput(args.transfo);
    trReader->SetInvertTransform(args.invert);
    trReader->Update();
    typename TransformType::Pointer transfo = trReader->GetOutputTransform();

    std::cout << "Image to transform is vector." << std::endl;

    typedef itk::ResampleImageFilter<ImageType, OuputType> ResampleFilterType;
    typename ResampleFilterType::Pointer vectorResampler = ResampleFilterType::New();
    typename itk::InterpolateImageFunction <ImageType>::Pointer interpolator;

    if(args.interpolation == "nearest")
        interpolator = itk::NearestNeighborInterpolateImageFunction<ImageType>::New();
    else if(args.interpolation == "linear")
        interpolator = itk::LinearInterpolateImageFunction<ImageType>::New();
    else
    {
        itk::ExceptionObject excp(__FILE__, __LINE__,
                                  "bspline and sinc interpolation not suported for vector images yet.",
                                  ITK_LOCATION);
        throw excp;
    }

    vectorResampler->SetTransform(transfo);
    vectorResampler->SetInterpolator(interpolator);

    typename ImageType::SizeType size;
    typename ImageType::PointType origin;
    typename ImageType::SpacingType spacing;
    typename ImageType::DirectionType direction;
    for(unsigned int i = 0; i < ImageType::ImageDimension; ++i)
    {
        size[i] = geometryImageIO->GetDimensions(i);
        origin[i] = geometryImageIO->GetOrigin(i);
        spacing[i] = geometryImageIO->GetSpacing(i);
        for(unsigned int j = 0; j < ImageType::ImageDimension; ++j)
            direction[i][j] = geometryImageIO->GetDirection(j)[i];
    }
    vectorResampler->SetSize(size);
    vectorResampler->SetOutputOrigin(origin);
    vectorResampler->SetOutputSpacing(spacing);
    vectorResampler->SetOutputDirection(direction);

    vectorResampler->SetInput(anima::readImage<ImageType>(args.input));
    vectorResampler->SetNumberOfThreads(args.pthread);
    vectorResampler->Update();

    anima::writeImage<OuputType>(args.output, vectorResampler->GetOutput());
}


template <class ImageType>
void
applyScalarTransfo(itk::ImageIOBase::Pointer geometryImageIO, const arguments &args)
{
    typedef itk::Image<float, ImageType::ImageDimension> OuputType;

    typedef anima::TransformSeriesReader <double, ImageType::ImageDimension> TransformSeriesReaderType;
    typedef typename TransformSeriesReaderType::OutputTransformType TransformType;
    TransformSeriesReaderType *trReader = new TransformSeriesReaderType;
    trReader->SetInput(args.transfo);
    trReader->SetInvertTransform(args.invert);
    trReader->Update();
    typename TransformType::Pointer transfo = trReader->GetOutputTransform();

    std::cout << "Image to transform is Scalar." << std::endl;
    typedef anima::ResampleImageFilter<ImageType, OuputType> ResampleFilterType;
    typename ResampleFilterType::Pointer scalarResampler = ResampleFilterType::New();
    typename itk::InterpolateImageFunction <ImageType>::Pointer interpolator;

    if(args.interpolation == "nearest")
        interpolator = itk::NearestNeighborInterpolateImageFunction<ImageType>::New();
    else if(args.interpolation == "linear")
        interpolator = anima::FasterLinearInterpolateImageFunction<ImageType>::New();
    else if(args.interpolation == "bspline")
        interpolator = itk::BSplineInterpolateImageFunction<ImageType>::New();
    else if(args.interpolation == "sinc")
    {
        const unsigned int WindowRadius = 4;
        typedef itk::Function::HammingWindowFunction<WindowRadius> WindowFunctionType;
        typedef itk::ConstantBoundaryCondition<ImageType> BoundaryConditionType;
        interpolator = itk::WindowedSincInterpolateImageFunction
                <ImageType, WindowRadius, WindowFunctionType, BoundaryConditionType, double >::New();
    }

    scalarResampler->SetTransform(transfo);
    scalarResampler->SetInterpolator(interpolator);

    typename ImageType::SizeType size;
    typename ImageType::PointType origin;
    typename ImageType::SpacingType spacing;
    typename ImageType::DirectionType direction;
    for(unsigned int i = 0; i < ImageType::ImageDimension; ++i)
    {
        size[i] = geometryImageIO->GetDimensions(i);
        origin[i] = geometryImageIO->GetOrigin(i);
        spacing[i] = geometryImageIO->GetSpacing(i);
        for(unsigned int j = 0; j < ImageType::ImageDimension; ++j)
            direction[i][j] = geometryImageIO->GetDirection(j)[i];
    }
    scalarResampler->SetSize(size);
    scalarResampler->SetOutputOrigin(origin);
    scalarResampler->SetOutputSpacing(spacing);
    scalarResampler->SetOutputDirection(direction);

    scalarResampler->SetInput(anima::readImage<ImageType>(args.input));
    scalarResampler->SetNumberOfThreads(args.pthread);
    scalarResampler->Update();

    anima::writeImage<OuputType>(args.output, scalarResampler->GetOutput());
}


template <class ComponentType, int Dimension>
void
checkIfComponentsAreVectors(itk::ImageIOBase::Pointer inputImageIO, itk::ImageIOBase::Pointer geometryImageIO, const arguments &args)
{
    if (inputImageIO->GetNumberOfComponents() > 1)
        applyVectorTransfo<itk::VectorImage<ComponentType, Dimension> >(geometryImageIO, args);
    else
        applyScalarTransfo<itk::Image<ComponentType, Dimension> >(geometryImageIO, args);
}

template <class ComponentType>
void
retrieveNbDimensions(itk::ImageIOBase::Pointer inputImageIO, itk::ImageIOBase::Pointer geometryImageIO,  const arguments &args)
{
    if(inputImageIO->GetNumberOfDimensions() != 3)
    {
        itk::ExceptionObject excp(__FILE__, __LINE__, "Number of type not supported.", ITK_LOCATION);
        throw excp;
    }
    else
        checkIfComponentsAreVectors<ComponentType, 3>(inputImageIO, geometryImageIO, args);
}


int main(int ac, const char** av)
{
    std::string descriptionMessage = "Resampler tool to apply args series of transformations to an image.\n"
                                     "Input transform is an XML file describing all transforms to apply.\n"
                                     "Such args file should look like this:\n"
                                     "<TransformationList>\n"
                                     "<Transformation>\n"
                                     "<Type>linear</Type> (it can be svf or dense too)\n"
                                     "<Path>FileName</Path>\n"
                                     "<Inversion>0</Inversion>\n"
                                     "</Transformation>\n"
                                     "...\n"
                                     "</TransformationList>\n"
                                     "Note that only geometries and input with the same number of dimensions are supported for now.\n"
                                     "The default interpolation methode is linear, it can be [nearest, linear, bspline, sinc]."
                                     "INRIA / IRISA - VisAGeS Team";

    TCLAP::CmdLine cmd(descriptionMessage, ' ',"1.0");

    TCLAP::ValueArg<std::string> inArg("i","input","Input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> trArg("t","trsf","Transformations XML list",true,"","transformations list",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","Output resampled image",true,"","output image",cmd);
    TCLAP::ValueArg<std::string> geomArg("g","geometry","Geometry image",true,"","geometry image",cmd);

    TCLAP::SwitchArg invertArg("I","invert","Invert the transformation series",cmd,false);
    TCLAP::ValueArg<std::string> interpolationArg("n",
                                                  "interpolation",
                                                  "interpolation method to use [nearest, linear, bspline, sinc]",
                                                  false,
                                                  "linear",
                                                  "interpolation method",
                                                  cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default : all cores)",
                                         false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(ac,av);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    // Find out the type of the image in file
    itk::ImageIOBase::Pointer inputImageIO = itk::ImageIOFactory::CreateImageIO(inArg.getValue().c_str(),
                                                                           itk::ImageIOFactory::ReadMode);
    itk::ImageIOBase::Pointer geometryImageIO = itk::ImageIOFactory::CreateImageIO(geomArg.getValue().c_str(),
                                                                           itk::ImageIOFactory::ReadMode);
    if(!inputImageIO)
    {
        std::cerr << "Itk could not find suitable IO factory for the input." << std::endl;
        return EXIT_FAILURE;
    }
    if(!geometryImageIO)
    {
        std::cerr << "Itk could not find suitable IO factory for the geometry image." << std::endl;
        return EXIT_FAILURE;
    }
    if(inputImageIO->GetNumberOfDimensions() != geometryImageIO->GetNumberOfDimensions())
    {
        std::cerr << "Input and geometry images have to have the same number of dimensions." << std::endl;
        return EXIT_FAILURE;
    }

    // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
    inputImageIO->SetFileName(inArg.getValue());
    inputImageIO->ReadImageInformation();
    geometryImageIO->SetFileName(geomArg.getValue());
    geometryImageIO->ReadImageInformation();

    arguments args;
    args.input = inArg.getValue();
    args.output = outArg.getValue();
    args.geometry = geomArg.getValue();
    args.transfo = trArg.getValue();
    args.invert = invertArg.getValue();
    args.pthread = nbpArg.getValue();
    args.interpolation = interpolationArg.getValue();

    bool badInterpolation = true;
    std::string interpolations[4] = {"nearest", "linear", "bspline", "sinc"};
    for(int i = 0; i < 4; ++i)
    {
        if(args.interpolation == interpolations[i])
        {
            badInterpolation = false;
            break;
        }
    }

    if(badInterpolation)
         std::cerr << "Interpolation method not suported, it must be one of [nearest, linear, bspline, sinc]." << std::endl;

    try
    {
        ANIMA_RETRIEVE_COMPONENT_TYPE(inputImageIO,
                                      retrieveNbDimensions,
                                      inputImageIO,
                                      geometryImageIO,
                                      args)
    }
    catch ( itk::ExceptionObject & err )
    {
        std::cerr << "Can't apply transformation, be sure to use valid arguments..." << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
