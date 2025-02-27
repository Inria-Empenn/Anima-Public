#include <tclap/CmdLine.h>
#include <itkWeightedCentroidKdTreeGenerator.h>
#include <itkKdTreeBasedKmeansEstimator.h>
#include <itkListSample.h>
#include <itkMeasurementVectorTraits.h>
#include <animaReadWriteFunctions.h>

void PerformLinearRegression(const std::vector<double> &x, const std::vector<double> &y,
                             double &slope, double &intercept)
{
    double correlation = 0.0;
    double meanX = 0.0;
    double meanY = 0.0;
    double varianceX = 0.0;
    unsigned int sizeOfArray = x.size();

    for (unsigned int arrayIndex = 0;arrayIndex < sizeOfArray;++arrayIndex)
    {
        meanX += x[arrayIndex];
        meanY += y[arrayIndex];
    }

    meanX /= sizeOfArray;
    meanY /= sizeOfArray;

    for (unsigned int arrayIndex = 0;arrayIndex < sizeOfArray;++arrayIndex)
    {
        correlation += (x[arrayIndex] - meanX) * (y[arrayIndex] - meanY);
        varianceX += (x[arrayIndex] - meanX) * (x[arrayIndex] - meanX);
    }

    slope = correlation / varianceX;
    intercept = meanY - slope * meanX;
}

template <typename PixelType>
void GetKMeansClassesMeans(itk::Image <PixelType, 3> *image, itk::Image <unsigned short, 3> *mask,
                           std::vector <double> &imageClassMeans)
{
    unsigned int numberOfClasses = 3;
    typedef itk::Image <PixelType, 3> ImageType;
    typedef itk::Image <unsigned short, 3> MaskType;

    itk::ImageRegionIterator <ImageType> imageItr(image, image->GetLargestPossibleRegion());
    itk::ImageRegionIterator <MaskType> maskItr(mask, mask->GetLargestPossibleRegion());

    unsigned int numPixels = 0;
    double imageMeanInMask = 0;
    while (!imageItr.IsAtEnd())
    {
        if (maskItr.Get() != 0)
        {
            imageMeanInMask += imageItr.Get();
            ++numPixels;
        }

        ++imageItr;
        ++maskItr;
    }

    imageMeanInMask /= numPixels;

    typedef itk::Statistics::MeasurementVectorPixelTraits <double>::MeasurementVectorType MeasurementVectorType;
    typedef itk::Statistics::ListSample <MeasurementVectorType> ListSampleType;
    typedef itk::Statistics::WeightedCentroidKdTreeGenerator <ListSampleType> TreeGeneratorType;

    ListSampleType::Pointer treeData = ListSampleType::New();
    treeData->Resize(numPixels);

    maskItr.GoToBegin();
    imageItr.GoToBegin();
    unsigned int pos = 0;
    while (!imageItr.IsAtEnd())
    {
        if (maskItr.Get() != 0)
        {
            treeData->SetMeasurementVector(pos, MeasurementVectorType(imageItr.Get()));
            ++pos;
        }

        ++imageItr;
        ++maskItr;
    }

    typedef TreeGeneratorType::KdTreeType TreeType;
    typedef itk::Statistics::KdTreeBasedKmeansEstimator <TreeType> EstimatorType;
    typedef EstimatorType::ParametersType ParametersType;

    TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
    treeGenerator->SetSample(treeData);
    treeGenerator->SetBucketSize(16);
    treeGenerator->Update();

    ParametersType initialMeans(numberOfClasses);
    initialMeans[0] = 0.3 * imageMeanInMask;
    initialMeans[1] = 0.7 * imageMeanInMask;
    initialMeans[2] = 0.9 * imageMeanInMask;

    EstimatorType::Pointer estimator = EstimatorType::New();
    estimator->SetParameters(initialMeans);

    estimator->SetKdTree(treeGenerator->GetOutput());
    estimator->SetMaximumIteration(200);
    estimator->SetCentroidPositionChangesThreshold(0.0);
    estimator->StartOptimization();

    ParametersType finalMeans = estimator->GetParameters();

    imageClassMeans.resize(numberOfClasses);
    for (unsigned int i = 0;i < numberOfClasses;++i)
        imageClassMeans[i] = finalMeans[i];
}

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd("Perform longitudinal intensity normalization. For more details on the method, see http://www.hal.inserm.fr/inserm-01074699/document", ' ', ANIMA_VERSION);

    TCLAP::ValueArg <std::string> referenceImageArg("r","ref","Reference image for intensity .",true,"","Reference image.",cmd);
    TCLAP::ValueArg <std::string> refMaskArg("R","ref-mask","Reference mask for intensity normalization.",true,"","Reference mask",cmd);
    TCLAP::ValueArg <std::string> movingImageArg("m","mov","Moving image for intensity normalization.",true,"","Moving image.", cmd);
    TCLAP::ValueArg <std::string> movingMaskArg("M","mov-mask","Moving mask for intensity normalization.",false,"","moving mask",cmd);
    TCLAP::ValueArg <std::string> outputArg("o", "output","Normalised output image",true,"","output image",cmd);

    TCLAP::ValueArg <unsigned int> numThreadsArg("T","num-threads","Number of cores to run on (default: all cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"Number of cores",cmd);

    try
    {
        cmd.parse(argc, argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef itk::Image <double, 3> ImageType;
    typedef itk::Image <unsigned short, 3> MaskType;

    ImageType::Pointer baselineImage = anima::readImage<ImageType>(referenceImageArg.getValue());
    ImageType::Pointer repeatImage = anima::readImage<ImageType>(movingImageArg.getValue());
    MaskType::Pointer refMask = anima::readImage<MaskType>(refMaskArg.getValue());

    MaskType::Pointer movingMask;
    if (movingMaskArg.getValue() != "")
        movingMask = anima::readImage<MaskType>(movingMaskArg.getValue());
    else
        movingMask = anima::readImage<MaskType>(refMaskArg.getValue());

    // Calculate the class means for each image
    unsigned int numberOfClasses = 3;
    std::vector <double> baselineFinalMeans(numberOfClasses);
    GetKMeansClassesMeans <double>(baselineImage, refMask, baselineFinalMeans);

    std::vector <double> repeatFinalMeans(numberOfClasses);
    GetKMeansClassesMeans <double>(repeatImage, movingMask, repeatFinalMeans);

    itk::ImageRegionIterator <ImageType> repeatItr(repeatImage, repeatImage->GetLargestPossibleRegion());
    itk::ImageRegionIterator <MaskType> movingMaskItr(movingMask, movingMask->GetLargestPossibleRegion());

    double slopeRepeatBaseline = 0.0;
    double interceptRepeatBaseline = 0.0;

    // Repeat (x) regress on baseline (y).
    PerformLinearRegression(repeatFinalMeans, baselineFinalMeans, slopeRepeatBaseline,
                            interceptRepeatBaseline);

    movingMaskItr.GoToBegin();
    for (repeatItr.GoToBegin(); !repeatItr.IsAtEnd(); ++repeatItr)
    {
        if (movingMaskItr.Get() != 0)
        {
            double repeatValue = repeatItr.Get();
            repeatItr.Set(std::abs(slopeRepeatBaseline * repeatValue + interceptRepeatBaseline));
        }

        ++movingMaskItr;
    }

    anima::writeImage<ImageType>(outputArg.getValue(), repeatImage);

    return EXIT_SUCCESS;
}
