#include <tclap/CmdLine.h>
#include <iostream>

#include <itkImage.h>
#include <itkResampleImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkWindowedSincInterpolateImageFunction.h>
#include <itkConstantBoundaryCondition.h>

#include <animaReadWriteFunctions.h>
#include <itkTranslationTransform.h>

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

    typedef itk::Image <float,3> ImageType;
    typedef itk::ResampleImageFilter <ImageType,ImageType> ResampleFilterType;

    ImageType::Pointer inputImage = anima::readImage <ImageType> (inArg.getValue());

    ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    resampler->SetInput(inputImage);

    ImageType::SizeType size = inputImage->GetLargestPossibleRegion().GetSize();
    ImageType::PointType origin = inputImage->GetOrigin();
    ImageType::SpacingType spacing = inputImage->GetSpacing();
    ImageType::DirectionType direction = inputImage->GetDirection();
    unsigned int dimension = ImageType::GetImageDimension();

    for (unsigned int i = 0;i < dimension;++i)
    {
        double outRes = (i == 0) * xArg.getValue() + (i == 1) * yArg.getValue() + (i == 2) * zArg.getValue();

        double spacingRatio = spacing[i] / outRes;
        double oldSize = size[i];
        size[i] = std::floor(size[i] * spacingRatio);
        double trueOutRes = oldSize * spacing[i] / size[i];
        origin[i] += 0.5 * (trueOutRes - spacing[i]);
        spacing[i] = trueOutRes;
    }

    resampler->SetSize(size);
    resampler->SetOutputOrigin(origin);
    resampler->SetOutputSpacing(spacing);
    resampler->SetOutputDirection(direction);

    itk::TranslationTransform <double, 3>::Pointer idTrsf = itk::TranslationTransform <double, 3>::New();
    idTrsf->SetIdentity();
    resampler->SetTransform(idTrsf);

    typename itk::InterpolateImageFunction <ImageType>::Pointer interpolator;

    if(interpolationArg.getValue() == "nearest")
        interpolator = itk::NearestNeighborInterpolateImageFunction<ImageType>::New();
    else if(interpolationArg.getValue() == "linear")
        interpolator = itk::LinearInterpolateImageFunction<ImageType>::New();
    else if(interpolationArg.getValue() == "bspline")
        interpolator = itk::BSplineInterpolateImageFunction<ImageType>::New();
    else if(interpolationArg.getValue() == "sinc")
    {
        const unsigned int WindowRadius = 4;
        typedef itk::Function::HammingWindowFunction<WindowRadius> WindowFunctionType;
        typedef itk::ConstantBoundaryCondition<ImageType> BoundaryConditionType;
        interpolator = itk::WindowedSincInterpolateImageFunction
                <ImageType, WindowRadius, WindowFunctionType, BoundaryConditionType, double >::New();
    }

    resampler->SetInterpolator(interpolator);
    resampler->SetNumberOfWorkUnits(nbpArg.getValue());
    resampler->Update();

    anima::writeImage <ImageType> (outArg.getValue(),resampler->GetOutput());

    return EXIT_SUCCESS;
}
