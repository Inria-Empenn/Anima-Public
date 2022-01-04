#include <tclap/CmdLine.h>
#include <iostream>

#include <animaResampleImageFilter.h>
#include <itkResampleImageFilter.h>

#include <itkImage.h>
#include <itkVectorImage.h>

#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkWindowedSincInterpolateImageFunction.h>
#include <itkConstantBoundaryCondition.h>

#include <animaRetrieveImageTypeMacros.h>
#include <animaReadWriteFunctions.h>
#include <itkTranslationTransform.h>

struct arguments
{
    unsigned int pthread;
    double xSpacing, ySpacing, zSpacing;
    std::string input, output, interpolation;
};

template <class ImageType>
void
changeScalarImageResolution4D(const arguments &args)
{
    typedef itk::Image <double, ImageType::ImageDimension> OutputType;
    const unsigned int InternalImageDimension = 3;
    typedef itk::Image <double, InternalImageDimension> InternalImageType;
    typename itk::InterpolateImageFunction <InternalImageType>::Pointer interpolator;

    if(args.interpolation == "nearest")
        interpolator = itk::NearestNeighborInterpolateImageFunction<InternalImageType>::New();
    else if(args.interpolation == "linear")
        interpolator = itk::LinearInterpolateImageFunction<InternalImageType>::New();
    else if(args.interpolation == "bspline")
        interpolator = itk::BSplineInterpolateImageFunction<InternalImageType>::New();
    else if(args.interpolation == "sinc")
    {
        const unsigned int WindowRadius = 4;
        typedef itk::Function::HammingWindowFunction<WindowRadius> WindowFunctionType;
        typedef itk::ConstantBoundaryCondition<InternalImageType> BoundaryConditionType;
        interpolator = itk::WindowedSincInterpolateImageFunction
                <InternalImageType, WindowRadius, WindowFunctionType, BoundaryConditionType, double >::New();
    }

    typename ImageType::Pointer inputImage = anima::readImage <ImageType> (args.input);

    typename OutputType::PointType origin;
    typename OutputType::SpacingType spacing;
    typename OutputType::DirectionType direction;
    direction.SetIdentity();

    typename InternalImageType::SizeType internalSize;
    typename InternalImageType::PointType internalOrigin;
    typename InternalImageType::SpacingType internalSpacing;
    typename InternalImageType::DirectionType internalDirection;

    typename OutputType::RegionType outputRegion;
    for (unsigned int i = 0;i < InternalImageDimension;++i)
    {
        double outRes = (i == 0) * args.xSpacing + (i == 1) * args.ySpacing + (i == 2) * args.zSpacing;

        double spacingRatio = inputImage->GetSpacing()[i] / outRes;
        double oldSize = inputImage->GetLargestPossibleRegion().GetSize()[i];
        unsigned int size = std::floor(oldSize * spacingRatio);
        double trueOutRes = oldSize * inputImage->GetSpacing()[i] / size;
        origin[i] = inputImage->GetOrigin()[i] + 0.5 * (trueOutRes - inputImage->GetSpacing()[i]);
        spacing[i] = trueOutRes;

        outputRegion.SetIndex(i,0);
        outputRegion.SetSize(i,size);
        for(unsigned int j = 0;j < InternalImageDimension;++j)
            direction(i,j) = inputImage->GetDirection()(i,j);

        internalSize[i] = size;
        internalOrigin[i] = origin[i];
        internalSpacing[i] = spacing[i];
        for(unsigned int j = 0;j < InternalImageDimension;++j)
            internalDirection(i,j) = direction(i,j);
    }

    for (unsigned int i = InternalImageDimension;i < ImageType::GetImageDimension();++i)
    {
        outputRegion.SetIndex(i,0);
        outputRegion.SetSize(i,inputImage->GetLargestPossibleRegion().GetSize()[i]);
        origin[i] = inputImage->GetOrigin()[i];
        spacing[i] = inputImage->GetSpacing()[i];
        direction(i,i) = inputImage->GetDirection()(i,i);
    }

    typename OutputType::Pointer outputImage = OutputType::New();
    outputImage->Initialize();
    outputImage->SetRegions(outputRegion);
    outputImage->SetOrigin(origin);
    outputImage->SetSpacing(spacing);
    outputImage->SetDirection(direction);
    outputImage->Allocate();

    unsigned int numImages = inputImage->GetLargestPossibleRegion().GetSize()[InternalImageDimension];

    itk::TranslationTransform <double, 3>::Pointer idTrsf = itk::TranslationTransform <double, 3>::New();
    idTrsf->SetIdentity();

    for (unsigned int i = 0;i < numImages;++i)
    {
        if (i == 0)
            std::cout << "Resampling sub-image " << i+1 << "/" << numImages << std::flush;
        else
            std::cout<<"\033[K\rResampling sub-image " << i+1 << "/" << numImages << std::flush;

        typedef anima::ResampleImageFilter<InternalImageType, InternalImageType> ResampleFilterType;
        typename ResampleFilterType::Pointer scalarResampler = ResampleFilterType::New();
        scalarResampler->SetTransform(idTrsf);
        scalarResampler->SetInterpolator(interpolator);

        scalarResampler->SetSize(internalSize);
        scalarResampler->SetOutputOrigin(internalOrigin);
        scalarResampler->SetOutputSpacing(internalSpacing);
        scalarResampler->SetOutputDirection(internalDirection);

        typedef itk::ExtractImageFilter <ImageType, InternalImageType> ExtractFilterType;
        typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
        extractFilter->SetInput(inputImage);
        extractFilter->SetDirectionCollapseToGuess();

        typename ImageType::RegionType extractRegion = inputImage->GetLargestPossibleRegion();
        extractRegion.SetIndex(InternalImageDimension,i);
        extractRegion.SetSize(InternalImageDimension,0);

        extractFilter->SetExtractionRegion(extractRegion);
        extractFilter->SetNumberOfWorkUnits(args.pthread);

        extractFilter->Update();

        scalarResampler->SetInput(extractFilter->GetOutput());
        scalarResampler->SetNumberOfWorkUnits(args.pthread);
        scalarResampler->Update();

        extractRegion.SetSize(InternalImageDimension,1);
        for (unsigned int j = 0;j < InternalImageDimension;++j)
        {
            extractRegion.SetIndex(j,scalarResampler->GetOutput()->GetLargestPossibleRegion().GetIndex()[j]);
            extractRegion.SetSize(j,scalarResampler->GetOutput()->GetLargestPossibleRegion().GetSize()[j]);
        }

        itk::ImageRegionIterator <InternalImageType> internalOutItr(scalarResampler->GetOutput(),
                                                                    scalarResampler->GetOutput()->GetLargestPossibleRegion());
        itk::ImageRegionIterator <OutputType> outItr(outputImage,extractRegion);

        while (!outItr.IsAtEnd())
        {
            outItr.Set(internalOutItr.Get());

            ++internalOutItr;
            ++outItr;
        }
    }

    std::cout << std::endl;

    anima::writeImage<OutputType>(args.output, outputImage);
}

template <class ImageType>
void
changeVectorImageResolution(const arguments &args)
{
    typedef itk::VectorImage <double,ImageType::ImageDimension> OutputType;
    typedef itk::ResampleImageFilter <ImageType,OutputType> ResampleFilterType;

    typename ImageType::Pointer inputImage = anima::readImage <ImageType> (args.input);

    typename ResampleFilterType::Pointer vectorResampler = ResampleFilterType::New();
    vectorResampler->SetInput(inputImage);

    typename ImageType::SizeType size = inputImage->GetLargestPossibleRegion().GetSize();
    typename ImageType::PointType origin = inputImage->GetOrigin();
    typename ImageType::SpacingType spacing = inputImage->GetSpacing();
    typename ImageType::DirectionType direction = inputImage->GetDirection();
    unsigned int dimension = ImageType::GetImageDimension();

    for (unsigned int i = 0;i < dimension;++i)
    {
        double outRes = (i == 0) * args.xSpacing + (i == 1) * args.ySpacing + (i == 2) * args.zSpacing;

        double spacingRatio = spacing[i] / outRes;
        double oldSize = size[i];
        size[i] = std::floor(oldSize * spacingRatio);
        double trueOutRes = oldSize * spacing[i] / size[i];
        origin[i] += 0.5 * (trueOutRes - spacing[i]);
        spacing[i] = trueOutRes;
    }

    vectorResampler->SetSize(size);
    vectorResampler->SetOutputOrigin(origin);
    vectorResampler->SetOutputSpacing(spacing);
    vectorResampler->SetOutputDirection(direction);

    itk::TranslationTransform <double, 3>::Pointer idTrsf = itk::TranslationTransform <double, 3>::New();
    idTrsf->SetIdentity();
    vectorResampler->SetTransform(idTrsf);

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

    vectorResampler->SetInterpolator(interpolator);
    vectorResampler->SetNumberOfWorkUnits(args.pthread);
    vectorResampler->Update();

    anima::writeImage<OutputType>(args.output, vectorResampler->GetOutput());
}

template <class ImageType>
void
changeScalarImageResolution(const arguments &args)
{
    typedef itk::Image <double,ImageType::ImageDimension> OutputType;
    typedef anima::ResampleImageFilter <ImageType,OutputType> ResampleFilterType;

    typename ImageType::Pointer inputImage = anima::readImage <ImageType> (args.input);

    typename ResampleFilterType::Pointer scalarResampler = ResampleFilterType::New();
    scalarResampler->SetInput(inputImage);

    typename ImageType::SizeType size = inputImage->GetLargestPossibleRegion().GetSize();
    typename ImageType::PointType origin = inputImage->GetOrigin();
    typename ImageType::SpacingType spacing = inputImage->GetSpacing();
    typename ImageType::DirectionType direction = inputImage->GetDirection();
    unsigned int dimension = ImageType::GetImageDimension();

    for (unsigned int i = 0;i < dimension;++i)
    {
        double outRes = (i == 0) * args.xSpacing + (i == 1) * args.ySpacing + (i == 2) * args.zSpacing;

        double spacingRatio = spacing[i] / outRes;
        double oldSize = size[i];
        size[i] = std::floor(size[i] * spacingRatio);
        double trueOutRes = oldSize * spacing[i] / size[i];
        origin[i] += 0.5 * (trueOutRes - spacing[i]);
        spacing[i] = trueOutRes;
    }

    scalarResampler->SetSize(size);
    scalarResampler->SetOutputOrigin(origin);
    scalarResampler->SetOutputSpacing(spacing);
    scalarResampler->SetOutputDirection(direction);

    itk::TranslationTransform <double, 3>::Pointer idTrsf = itk::TranslationTransform <double, 3>::New();
    idTrsf->SetIdentity();
    scalarResampler->SetTransform(idTrsf);

    typename itk::InterpolateImageFunction <ImageType>::Pointer interpolator;

    if(args.interpolation == "nearest")
        interpolator = itk::NearestNeighborInterpolateImageFunction<ImageType>::New();
    else if(args.interpolation == "linear")
        interpolator = itk::LinearInterpolateImageFunction<ImageType>::New();
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

    scalarResampler->SetInterpolator(interpolator);
    scalarResampler->SetNumberOfWorkUnits(args.pthread);
    scalarResampler->Update();

    anima::writeImage<OutputType>(args.output, scalarResampler->GetOutput());
}

template <class ComponentType, int Dimension>
void
checkIfComponentsAreVectors(itk::ImageIOBase::Pointer inputImageIO, const arguments &args)
{
    if (inputImageIO->GetNumberOfComponents() > 1)
    {
        if (Dimension > 3)
            throw itk::ExceptionObject (__FILE__, __LINE__, "Number of dimensions not supported for vector image resampling", ITK_LOCATION);

        changeVectorImageResolution < itk::VectorImage<ComponentType, 3> > (args);
    }
    else
    {
        if (Dimension < 4)
            changeScalarImageResolution < itk::Image<ComponentType, 3> > (args);
        else
            changeScalarImageResolution4D < itk::Image<ComponentType, 4> > (args);
    }
}

template <class ComponentType>
void
retrieveNbDimensions(itk::ImageIOBase::Pointer inputImageIO, const arguments &args)
{
    if (inputImageIO->GetNumberOfDimensions() > 4)
        throw itk::ExceptionObject(__FILE__, __LINE__, "Number of dimensions not supported.", ITK_LOCATION);

    if (inputImageIO->GetNumberOfDimensions() > 3)
        checkIfComponentsAreVectors<ComponentType, 4> (inputImageIO, args);
    else
        checkIfComponentsAreVectors<ComponentType, 3> (inputImageIO, args);
}

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","input","Input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","Output resampled image",true,"","output image",cmd);

    TCLAP::ValueArg<double> xArg("x","x-res","Output voxel spacing on X direction (default: 1.0)",false,1.0,"voxel spacing on X",cmd);
    TCLAP::ValueArg<double> yArg("y","y-res","Output voxel spacing on Y direction (default: 1.0)",false,1.0,"voxel spacing on Y",cmd);
    TCLAP::ValueArg<double> zArg("z","z-res","Output voxel spacing on Z direction (default: 1.0)",false,1.0,"voxel spacing on Z",cmd);

    TCLAP::ValueArg<std::string> interpolationArg("n","interpolation","interpolation method to use [nearest, linear, bspline, sinc]",
                                                  false,"linear","interpolation method",cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default : all cores)",
                                         false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    // Find out the type of the image in file
    itk::ImageIOBase::Pointer inputImageIO = itk::ImageIOFactory::CreateImageIO(inArg.getValue().c_str(),
                                                                                itk::IOFileModeEnum::ReadMode);
    if(!inputImageIO)
    {
        std::cerr << "Itk could not find suitable IO factory for the input." << std::endl;
        return EXIT_FAILURE;
    }

    // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
    inputImageIO->SetFileName(inArg.getValue());
    inputImageIO->ReadImageInformation();

    arguments args;
    args.input = inArg.getValue();
    args.output = outArg.getValue();
    args.pthread = nbpArg.getValue();
    args.interpolation = interpolationArg.getValue();
    args.xSpacing = xArg.getValue();
    args.ySpacing = yArg.getValue();
    args.zSpacing = zArg.getValue();

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

    if (badInterpolation)
    {
        std::cerr << "Interpolation method not supported, it must be one of [nearest, linear, bspline, sinc]." << std::endl;
        return EXIT_FAILURE;
    }

    try
    {
        ANIMA_RETRIEVE_COMPONENT_TYPE(inputImageIO,
                                      retrieveNbDimensions,
                                      inputImageIO,
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
