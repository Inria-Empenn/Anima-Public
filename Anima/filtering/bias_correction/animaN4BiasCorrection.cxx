#include <tclap/CmdLine.h>
#include <itkN4BiasFieldCorrectionImageFilter.h>
#include <animaReadWriteFunctions.h>
#include <itkConstantPadImageFilter.h>
#include <itkOtsuThresholdImageFilter.h>
#include <itkShrinkImageFilter.h>
#include <itkCommand.h>
#include <itkTimeProbe.h>

//Update progression of the process
void eventCallback(itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
    std::cout << "\033[K\rProgression: " << (int)(processObject->GetProgress() * 100) << "%" << std::flush;
}

template<typename T>
T convertStringToNumber(std::string const &pi_rsValue)
{
    T numRes;

    std::istringstream(pi_rsValue)>> numRes;

    return numRes;
}

template<typename T>
std::vector<T> ConvertVector(std::string &pi_rsVector)
{
    std::vector<T>    oValuesVectorRes;

    size_t iSeparatorPosition = pi_rsVector.find('x', 0);

    if (iSeparatorPosition == std::string::npos)
    {
        oValuesVectorRes.push_back(convertStringToNumber<T>(pi_rsVector));
    }
    else
    {
        std::string sChunk = pi_rsVector.substr(0, iSeparatorPosition);

        oValuesVectorRes.push_back(convertStringToNumber<T>(sChunk));
        while (iSeparatorPosition != std::string::npos)
        {
            std::string::size_type crossposfrom = iSeparatorPosition;
            iSeparatorPosition = pi_rsVector.find('x', crossposfrom + 1);
            if (iSeparatorPosition == std::string::npos)
            {
                sChunk = pi_rsVector.substr(crossposfrom + 1, pi_rsVector.length());
            }
            else
            {
                sChunk = pi_rsVector.substr(crossposfrom + 1, iSeparatorPosition);
            }
            oValuesVectorRes.push_back(convertStringToNumber<unsigned int>(sChunk));
        }
    }

    return oValuesVectorRes;
}

void checkIterationArg(std::string &sIterations)
{
    char const *pchTemp = sIterations.c_str();
    for (int i = 0; i < sIterations.length(); ++i)
    {
        if (!((pchTemp[i] >= '0' && pchTemp[i] <= '9') || pchTemp[i] == 'x'))
        {
            std::cerr << "error: " << " for argument iterations does not respect format strictly positive number separated by 'x' like 100x50x50..." << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}

std::vector<unsigned int> extractMaxNumberOfIterationsVector(std::string pi_rsIterations)
{
    checkIterationArg(pi_rsIterations);
    return ConvertVector<unsigned int>(pi_rsIterations);
}

int main(int argc, char *argv[])
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ', ANIMA_VERSION);
    TCLAP::ValueArg<std::string> oArgInputImg("i", "input", "Input image.", true, "", "input Name", cmd);
    TCLAP::ValueArg<std::string> oArgOutputName("o", "output", "Name for output file", true, "", "output Name", cmd);
    TCLAP::ValueArg<std::string> oArgIterationsString("I", "iterations", "Table of number of iterations, default=50x40x30", false, "50x40x30", "iterations table", cmd);
    TCLAP::ValueArg<int> oArgShrinkFactors("S", "shrinkFactor", "Shrink factor, default=4", false, 4, "shrink Factor", cmd);
    TCLAP::ValueArg<float> oArgWienerFilterNoise("W", "wiener", "Wiener Filter Noise, default=0.01", false, 0.01, "wiener noise", cmd);
    TCLAP::ValueArg<float> oArgbfFWHM("B", "bfFWHM", "Bias field Full Width at Half Maximum, default=0.15", false, 0.15, "Bias field Full Width at Half Maximum", cmd);
    TCLAP::ValueArg<float> oArgConvergenceThreshold("T", "threshold", "Convergence Threshold, default=0.0001", false, 0.0001, "threshold", cmd);
    TCLAP::ValueArg<int> oArgSplineOrder("O", "splineOrder", "BSpline Order, default=3", false, 3, "spline Order", cmd);
    TCLAP::ValueArg<float> oArgSplineDistance("D", "splineDistance", "B-Spline distance, default=0.0", false, 0.0, "spline Distance", cmd);
    TCLAP::ValueArg<std::string> oArgInitialMeshResolutionString("G", "splineGrid", "B-Spline grid resolution. It is ignored if splineDistance>0 or if dimention of it <>3, default=1x1x1", false, "1x1x1", "spline Grid", cmd);
    TCLAP::ValueArg<unsigned int> oArgNbP("p", "numberofthreads", "Number of threads to run on (default: all cores)", false, itk::MultiThreader::GetGlobalDefaultNumberOfThreads(), "number of threads", cmd);

    try
    {
        cmd.parse( argc, argv );
    }
    catch (TCLAP::ArgException & e)
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
        exit(EXIT_FAILURE);
    }

    typedef itk::Image<float, 3 > ImageType;
    typedef itk::Image<unsigned char, 3> MaskImageType;
    typedef itk::N4BiasFieldCorrectionImageFilter<ImageType, MaskImageType, ImageType> BiasFilter;
    typedef itk::ConstantPadImageFilter<ImageType, ImageType> PadderType;
    typedef itk::ConstantPadImageFilter<MaskImageType, MaskImageType> MaskPadderType;
    typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkerType;
    typedef itk::ShrinkImageFilter<MaskImageType, MaskImageType> MaskShrinkerType;
    typedef itk::BSplineControlPointImageFilter< BiasFilter::BiasFieldControlPointLatticeType, BiasFilter::ScalarImageType> BSplinerType;
    typedef itk::ExpImageFilter<ImageType, ImageType> ExpFilterType;
    typedef itk::DivideImageFilter<ImageType, ImageType, ImageType> DividerType;
    typedef itk::ExtractImageFilter<ImageType, ImageType> CropperType;

    std::vector<unsigned int> oMaxNumbersIterationsVector = extractMaxNumberOfIterationsVector(oArgIterationsString.getValue());
    std::vector<float> oInitialMeshResolutionVect = ConvertVector<float>(oArgInitialMeshResolutionString.getValue());
    float fSplineDistance = oArgSplineDistance.getValue();


    /********************************************************************************/
    /***************************** PREPARING STARTING *******************************/
    /********************************************************************************/

    /*** 0 ******************* Create filter and accessories ******************/
    BiasFilter::Pointer filter = BiasFilter::New();
    BiasFilter::ArrayType oNumberOfControlPointsArray;

    /*** 1 ******************* Read input image *******************************/
    ImageType::Pointer image = anima::readImage<ImageType>( oArgInputImg.getValue());

    /*** 2 ******************* Creating Otsu mask *****************************/
    std::cout << "Creating Otsu mask." << std::endl;
    itk::TimeProbe timer;
    timer.Start();
    MaskImageType::Pointer maskImage = ITK_NULLPTR;
    typedef itk::OtsuThresholdImageFilter<ImageType, MaskImageType> ThresholderType;
    ThresholderType::Pointer otsu = ThresholderType::New();
    otsu->SetInput(image);
    otsu->SetNumberOfHistogramBins(200);
    otsu->SetInsideValue(0);
    otsu->SetOutsideValue(1);

    otsu->SetNumberOfThreads(oArgNbP.getValue());
    otsu->Update();
    maskImage = otsu->GetOutput();

    /*** 3A *************** Set Maximum number of Iterations for the filter ***/
    BiasFilter::VariableSizeArrayType itkTabMaximumIterations;
    itkTabMaximumIterations.SetSize(oMaxNumbersIterationsVector.size());
    for (int i=0; i<oMaxNumbersIterationsVector.size(); ++i)
    {
        itkTabMaximumIterations[i] = oMaxNumbersIterationsVector[i];
    }
    filter->SetMaximumNumberOfIterations(itkTabMaximumIterations);

    /*** 3B *************** Set Fitting Levels for the filter *****************/
    BiasFilter::ArrayType oFittingLevelsTab;
    oFittingLevelsTab.Fill(oMaxNumbersIterationsVector.size());
    filter->SetNumberOfFittingLevels(oFittingLevelsTab);

    /*** 4 ******************* Save image's index, size, origine **************/
    ImageType::IndexType oImageIndex = image->GetLargestPossibleRegion().GetIndex();
    ImageType::SizeType oImageSize = image->GetLargestPossibleRegion().GetSize();
    ImageType::PointType newOrigin = image->GetOrigin();

    if (fSplineDistance>0)
    {
        /*** 5 ******************* Compute number of control points  **************/
        itk::SizeValueType lowerBound[3];
        itk::SizeValueType upperBound[3];

        for (unsigned int i = 0; i < 3; i++)
        {
            float domain = static_cast<float>(image->GetLargestPossibleRegion().GetSize()[i] - 1) * image->GetSpacing()[i];
            unsigned int numberOfSpans = static_cast<unsigned int>(std::ceil(domain / fSplineDistance));
            unsigned long extraPadding = static_cast<unsigned long>((numberOfSpans * fSplineDistance - domain) / image->GetSpacing()[i] + 0.5);
            lowerBound[i] = static_cast<unsigned long>(0.5 * extraPadding);
            upperBound[i] = extraPadding - lowerBound[i];
            newOrigin[i] -= (static_cast<float>(lowerBound[i]) * image->GetSpacing()[i]);
            oNumberOfControlPointsArray[i] = numberOfSpans + filter->GetSplineOrder();
        }

        /*** 6 ******************* Padder  ****************************************/
        PadderType::Pointer imagePadder = PadderType::New();
        imagePadder->SetInput(image);
        imagePadder->SetPadLowerBound(lowerBound);
        imagePadder->SetPadUpperBound(upperBound);
        imagePadder->SetConstant(0);
        imagePadder->SetNumberOfThreads(oArgNbP.getValue());
        imagePadder->Update();

        image = imagePadder->GetOutput();

        /*** 7 ******************** Handle the mask image *************************/
        MaskPadderType::Pointer maskPadder = MaskPadderType::New();
        maskPadder->SetInput(maskImage);
        maskPadder->SetPadLowerBound(lowerBound);
        maskPadder->SetPadUpperBound(upperBound);
        maskPadder->SetConstant(0);

        maskPadder->SetNumberOfThreads(oArgNbP.getValue());
        maskPadder->Update();

        maskImage = maskPadder->GetOutput();

        /*** 8 ******************** SetNumber Of Control Points *******************/
        filter->SetNumberOfControlPoints(oNumberOfControlPointsArray);
    }
    else if(oInitialMeshResolutionVect.size() == 3)
    {
        /*** 9 ******************** SetNumber Of Control Points alternative *******/
        for (unsigned i = 0; i < 3; i++)
        {
            oNumberOfControlPointsArray[i] = static_cast<unsigned int>(oInitialMeshResolutionVect[i]) + filter->GetSplineOrder();
        }
        filter->SetNumberOfControlPoints(oNumberOfControlPointsArray);
    }
    else
    {
        std::cout << "No BSpline distance and Mesh Resolution is ignored because not 3 dimensions" << std::endl;
    }

    /*** 10 ******************* Shrinker image ********************************/
    ShrinkerType::Pointer imageShrinker = ShrinkerType::New();
    imageShrinker->SetInput(image);
    imageShrinker->SetShrinkFactors(1);

    /*** 11 ******************* Shrinker mask *********************************/
    MaskShrinkerType::Pointer maskShrinker = MaskShrinkerType::New();
    maskShrinker->SetInput(maskImage);
    maskShrinker->SetShrinkFactors(1);

    /*** 12 ******************* Shrink mask and image *************************/
    imageShrinker->SetShrinkFactors(oArgShrinkFactors.getValue());
    maskShrinker->SetShrinkFactors(oArgShrinkFactors.getValue());
    imageShrinker->SetNumberOfThreads(oArgNbP.getValue());
    maskShrinker->SetNumberOfThreads(oArgNbP.getValue());
    imageShrinker->Update();
    maskShrinker->Update();

    /*** 13 ******************* Filter setings ********************************/
    filter->SetSplineOrder(oArgSplineOrder.getValue());
    filter->SetWienerFilterNoise(oArgWienerFilterNoise.getValue());
    filter->SetBiasFieldFullWidthAtHalfMaximum(oArgbfFWHM.getValue());
    filter->SetConvergenceThreshold(oArgConvergenceThreshold.getValue());
    filter->SetInput(imageShrinker->GetOutput());
    filter->SetMaskImage(maskShrinker->GetOutput());

    /*** 14 ******************* Apply filter **********************************/
    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);
    filter->AddObserver(itk::ProgressEvent(), callback);
    try
    {
        filter->SetNumberOfThreads(oArgNbP.getValue());
        filter->Update();
    }
    catch (itk::ExceptionObject & err)
    {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        exit(EXIT_FAILURE);
    }


    /**
   * Reconstruct the bias field at full image resolution.  Divide
   * the original input image by the bias field to get the final
   * corrected image.
   */
    BSplinerType::Pointer bspliner = BSplinerType::New();
    bspliner->SetInput(filter->GetLogBiasFieldControlPointLattice());
    bspliner->SetSplineOrder(filter->GetSplineOrder());
    bspliner->SetSize(image->GetLargestPossibleRegion().GetSize());
    bspliner->SetOrigin(newOrigin);
    bspliner->SetDirection(image->GetDirection());
    bspliner->SetSpacing(image->GetSpacing());
    bspliner->SetNumberOfThreads(oArgNbP.getValue());
    bspliner->Update();

    ImageType::Pointer logField = ImageType::New();
    logField->SetOrigin(image->GetOrigin());
    logField->SetSpacing(image->GetSpacing());
    logField->SetRegions(image->GetLargestPossibleRegion());
    logField->SetDirection(image->GetDirection());
    logField->Allocate();

    itk::ImageRegionIterator<BiasFilter::ScalarImageType> IB(bspliner->GetOutput(), bspliner->GetOutput()->GetLargestPossibleRegion());

    itk::ImageRegionIterator<ImageType> IF(logField, logField->GetLargestPossibleRegion());

    for (IB.GoToBegin(), IF.GoToBegin(); !IB.IsAtEnd(); ++IB, ++IF)
    {
        IF.Set(IB.Get()[0]);
    }

    ExpFilterType::Pointer expFilter = ExpFilterType::New();
    expFilter->SetInput(logField);
    expFilter->SetNumberOfThreads(oArgNbP.getValue());
    expFilter->Update();

    DividerType::Pointer divider = DividerType::New();
    divider->SetInput1(image);
    divider->SetInput2(expFilter->GetOutput());
    divider->SetNumberOfThreads(oArgNbP.getValue());
    divider->Update();

    ImageType::RegionType inputRegion;
    inputRegion.SetIndex(oImageIndex);
    inputRegion.SetSize(oImageSize);

    CropperType::Pointer cropper = CropperType::New();
    cropper->SetInput(divider->GetOutput());
    cropper->SetExtractionRegion(inputRegion);
    cropper->SetDirectionCollapseToSubmatrix();
    cropper->SetNumberOfThreads(oArgNbP.getValue());
    cropper->Update();

    timer.Stop();
    std::cout << "\nComputation time : " << timer.GetTotal() << std::endl;

    /********************** Write output image *************************/
    try
    {
        anima::writeImage<ImageType>(oArgOutputName.getValue(), cropper->GetOutput());
    }
    catch (itk::ExceptionObject & err)
    {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}

