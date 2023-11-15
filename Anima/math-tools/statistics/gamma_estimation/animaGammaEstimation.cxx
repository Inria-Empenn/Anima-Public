#include <animaGammaDistribution.h>
#include <animaReadWriteFunctions.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ', ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg(
        "i", "input-files",
        "Input image list in text file",
        true, "", "input image list", cmd);
    TCLAP::ValueArg<std::string> kappaArg(
        "k", "kappa-file",
        "Output kappa image",
        true, "", "output kappa image", cmd);
    TCLAP::ValueArg<std::string> thetaArg(
        "t", "theta-file",
        "Output theta image",
        true, "", "output theta image", cmd);
    TCLAP::ValueArg<std::string> maskArg(
        "m", "mask-files",
        "Input masks list in text file (mask images should contain only zeros or ones)",
        false, "", "input masks list", cmd);
    TCLAP::ValueArg<std::string> meanArg(
        "", "mean-file",
        "A 3D scalar image storing the mean of the Gamma prior.",
        false, "", "output mean image", cmd);
    TCLAP::ValueArg<std::string> varArg(
        "", "variance-file",
        "A 3D scalar image storing the variance of the Gamma prior.",
        false, "", "output variance image", cmd);
    TCLAP::ValueArg<std::string> methodArg(
        "", "method",
        "Estimation method. Choices are: mle, biased-closed-form or unbiased-closed-form",
        true, "mle", "estimation method", cmd);

    try
    {
        cmd.parse(argc, argv);
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    using InputImageType = itk::Image<double, 3>;
    using MaskImageType = itk::Image<unsigned int, 3>;
    using OutputImageType = itk::Image<double, 3>;

    // Read input sample
    std::ifstream imageIn(inArg.getValue());
    char imageN[2048];
    std::vector<InputImageType::Pointer> inputImages;
    while (!imageIn.eof())
    {
        imageIn.getline(imageN, 2048);
        if (strcmp(imageN, "") == 0)
            continue;
        inputImages.push_back(anima::readImage<InputImageType>(imageN));
    }
    imageIn.close();
    unsigned int nbImages = inputImages.size();

    std::vector<MaskImageType::Pointer> maskImages;
    if (maskArg.getValue() != "")
    {
        // Read input masks
        std::ifstream masksIn(maskArg.getValue());
        char maskN[2048];
        while (!masksIn.eof())
        {
            masksIn.getline(maskN, 2048);
            if (strcmp(maskN, "") == 0)
                continue;
            maskImages.push_back(anima::readImage<MaskImageType>(maskN));
        }
        masksIn.close();

        if (maskImages.size() != nbImages)
        {
            std::cerr << "The number of mask images should match the number of input images." << std::endl;
            return EXIT_FAILURE;
        }
    }

    using InputImageIteratorType = itk::ImageRegionConstIterator<InputImageType>;
    using MaskImageIteratorType = itk::ImageRegionConstIterator<MaskImageType>;
    using OutputImageIteratorType = itk::ImageRegionIterator<OutputImageType>;
    std::vector<InputImageIteratorType> inItrs(nbImages);
    std::vector<MaskImageIteratorType> maskItrs(nbImages);
    for (unsigned int i = 0; i < nbImages; ++i)
    {
        inItrs[i] = InputImageIteratorType(inputImages[i], inputImages[i]->GetLargestPossibleRegion());
        if (maskArg.getValue() != "")
            maskItrs[i] = MaskImageIteratorType(maskImages[i], maskImages[i]->GetLargestPossibleRegion());
    }

    // Initialize kappa image
    OutputImageType::Pointer kappaImage = OutputImageType::New();
    kappaImage->SetRegions(inputImages[0]->GetLargestPossibleRegion());
    kappaImage->CopyInformation(inputImages[0]);
    kappaImage->Allocate();
    kappaImage->FillBuffer(0.0);
    OutputImageIteratorType kappaItr(kappaImage, kappaImage->GetLargestPossibleRegion());

    // Initialize theta image
    OutputImageType::Pointer thetaImage = OutputImageType::New();
    thetaImage->SetRegions(inputImages[0]->GetLargestPossibleRegion());
    thetaImage->CopyInformation(inputImages[0]);
    thetaImage->Allocate();
    thetaImage->FillBuffer(0.0);
    OutputImageIteratorType thetaItr(thetaImage, thetaImage->GetLargestPossibleRegion());

    OutputImageType::Pointer meanImage;
    OutputImageIteratorType meanItr;
    if (meanArg.getValue() != "")
    {
        meanImage = OutputImageType::New();
        meanImage->SetRegions(inputImages[0]->GetLargestPossibleRegion());
        meanImage->CopyInformation(inputImages[0]);
        meanImage->Allocate();
        meanImage->FillBuffer(0.0);
        meanItr = OutputImageIteratorType(meanImage, meanImage->GetLargestPossibleRegion());
    }

    OutputImageType::Pointer varImage;
    OutputImageIteratorType varItr;
    if (varArg.getValue() != "")
    {
        varImage = OutputImageType::New();
        varImage->SetRegions(inputImages[0]->GetLargestPossibleRegion());
        varImage->CopyInformation(inputImages[0]);
        varImage->Allocate();
        varImage->FillBuffer(0.0);
        varItr = OutputImageIteratorType(varImage, varImage->GetLargestPossibleRegion());
    }

    std::vector<double> inputValues;
    anima::GammaDistribution gammaDistribution;
    double epsValue = std::sqrt(std::numeric_limits<double>::epsilon());
    std::vector<bool> usefulValues(nbImages, false);

    while (!kappaItr.IsAtEnd())
    {
        unsigned int nbUsedImages = 0;
        std::fill(usefulValues.begin(), usefulValues.end(), false);
        if (maskArg.getValue() != "")
        {
            for (unsigned int i = 0; i < nbImages; ++i)
            {
                if (maskItrs[i].Get() && inItrs[i].Get() > epsValue)
                {
                    usefulValues[i] = true;
                    nbUsedImages++;
                }
            }
        }
        else
        {
            for (unsigned int i = 0; i < nbImages; ++i)
            {
                if (inItrs[i].Get() > epsValue)
                {
                    usefulValues[i] = true;
                    nbUsedImages++;
                }
            }
        }

        if (nbUsedImages < 2)
        {
            for (unsigned int i = 0; i < nbImages; ++i)
            {
                ++inItrs[i];
                if (maskArg.getValue() != "")
                    ++maskItrs[i];
            }
            ++kappaItr;
            ++thetaItr;
            if (meanImage)
                ++meanItr;
            if (varImage)
                ++varItr;
            continue;
        }

        inputValues.resize(nbUsedImages);
        double meanValue = 0.0;
        unsigned int pos = 0;
        for (unsigned int i = 0; i < nbImages; ++i)
        {
            if (usefulValues[i])
            {
                double tmpValue = inItrs[i].Get();
                inputValues[pos] = tmpValue;
                meanValue += tmpValue;
                ++pos;
            }
        }

        meanValue /= static_cast<double>(nbUsedImages);
        double varValue = 0.0;
        for (unsigned int i = 0; i < nbUsedImages; ++i)
            varValue += (inputValues[i] - meanValue) * (inputValues[i] - meanValue);
        varValue /= (nbUsedImages - 1.0);

        if (varValue < epsValue)
        {
            for (unsigned int i = 0; i < nbImages; ++i)
            {
                ++inItrs[i];
                if (maskArg.getValue() != "")
                    ++maskItrs[i];
            }
            ++kappaItr;
            ++thetaItr;
            if (meanImage)
                ++meanItr;
            if (varImage)
                ++varItr;
            continue;
        }

        gammaDistribution.Fit(inputValues, methodArg.getValue());

        kappaItr.Set(gammaDistribution.GetShapeParameter());
        thetaItr.Set(gammaDistribution.GetScaleParameter());

        if (meanImage)
            meanItr.Set(gammaDistribution.GetMean());
        if (varImage)
            varItr.Set(gammaDistribution.GetVariance());

        for (unsigned int i = 0; i < nbImages; ++i)
        {
            ++inItrs[i];
            if (maskArg.getValue() != "")
                ++maskItrs[i];
        }
        ++kappaItr;
        ++thetaItr;
        if (meanImage)
            ++meanItr;
        if (varImage)
            ++varItr;
    }

    anima::writeImage<OutputImageType>(kappaArg.getValue(), kappaImage);
    anima::writeImage<OutputImageType>(thetaArg.getValue(), thetaImage);

    if (meanImage)
        anima::writeImage<OutputImageType>(meanArg.getValue(), meanImage);

    if (varImage)
        anima::writeImage<OutputImageType>(varArg.getValue(), varImage);

    return EXIT_SUCCESS;
}