#include <animaReadWriteFunctions.h>

#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkMultiThreader.h>
#include <vnl/vnl_matrix.h>

#include <tclap/CmdLine.h>
#include <fstream>
#include <algorithm>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ', ANIMA_VERSION);

    TCLAP::ValueArg <std::string> refArg("r", "ref", "Reference binary segmentation", true, "", "reference segmentation", cmd);
    TCLAP::ValueArg <std::string> testArg("t", "test", "Test binary segmentation", true, "", "test segmentation", cmd);
    TCLAP::ValueArg <std::string> outArg("o", "out", "output image with labeled lesions", true, "", "output labeled lesions", cmd);

    TCLAP::ValueArg <double> alphaArg("a", "alpha", "Alpha threshold (in between 0 and 1, default: 0)", false, 0, "alpha threshold", cmd);
    TCLAP::ValueArg <double> gammaArg("g", "gamma", "Gamma threshold (*100 %, default: 0.5)", false, 0.5, "gamma threshold", cmd);
    TCLAP::ValueArg <double> betaArg("b", "beta", "Beta threshold (in between 0 and 1, default: 0.05)", false, 0.05, "beta threshold", cmd);

    TCLAP::ValueArg <unsigned int> minVolumeArg("m", "min-vol", "Minimal volume for the component to be considered (default: 3 mm3)", false, 3, "minimal volume", cmd);
    TCLAP::ValueArg <unsigned int> nbpArg("T","numberofthreads","Number of threads to run on (default : all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    TCLAP::SwitchArg fullConnectArg("F","full-connect","Use 26-connectivity instead of 6-connectivity",cmd,false);

    try
    {
        cmd.parse(argc, argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef itk::Image <unsigned short, 3> ImageType;
    typedef itk::ImageRegionIterator <ImageType> ImageIteratorType;

    ImageType::Pointer refSegmentation = anima::readImage <ImageType> (refArg.getValue());
    ImageType::Pointer testSegmentation = anima::readImage <ImageType> (testArg.getValue());

    typedef itk::ConnectedComponentImageFilter <ImageType, ImageType> CCFilterType;
    typedef itk::RelabelComponentImageFilter <ImageType, ImageType> RelabelComponentFilterType;

    CCFilterType::Pointer refCCFilter = CCFilterType::New();
    refCCFilter->SetInput(refSegmentation);
    refCCFilter->SetFullyConnected(fullConnectArg.isSet());
    refCCFilter->SetNumberOfThreads(nbpArg.getValue());
    refCCFilter->Update();

    ImageType::SpacingType spacing = refSegmentation->GetSpacing();
    ImageType::SpacingValueType spacingTot = spacing[0];
    for (unsigned int i = 1;i < 3;++i)
        spacingTot *= spacing[i];

    // Compute minsize in voxels
    unsigned int minSizeInVoxel = static_cast <unsigned int> (std::ceil(minVolumeArg.getValue() / spacingTot));

    // Remove too small reference objects
    RelabelComponentFilterType::Pointer relabelRefFilter = RelabelComponentFilterType::New();
    relabelRefFilter->SetInput(refCCFilter->GetOutput());
    relabelRefFilter->SetMinimumObjectSize(minSizeInVoxel);
    relabelRefFilter->SetNumberOfThreads(nbpArg.getValue());
    relabelRefFilter->Update();

    // Reference segmentation is now labeled per connected objects
    refSegmentation = relabelRefFilter->GetOutput();
    refSegmentation->DisconnectPipeline();

    CCFilterType::Pointer testCCFilter = CCFilterType::New();
    testCCFilter->SetInput(testSegmentation);
    testCCFilter->SetFullyConnected(fullConnectArg.isSet());
    testCCFilter->SetNumberOfThreads(nbpArg.getValue());
    testCCFilter->Update();

    // Remove too small test objects
    RelabelComponentFilterType::Pointer relabelTestFilter = RelabelComponentFilterType::New();
    relabelTestFilter->SetInput(testCCFilter->GetOutput());
    relabelTestFilter->SetMinimumObjectSize(minSizeInVoxel);
    relabelTestFilter->SetNumberOfThreads(nbpArg.getValue());
    relabelTestFilter->Update();

    // Test segmentation is now labeled per connected objects
    testSegmentation = relabelTestFilter->GetOutput();
    testSegmentation->DisconnectPipeline();

    ImageIteratorType testItr(testSegmentation, testSegmentation->GetLargestPossibleRegion());

    unsigned short maxTestLabel = 0;
    while (!testItr.IsAtEnd())
    {
        if (testItr.Get() > maxTestLabel)
            maxTestLabel = testItr.Get();

        ++testItr;
    }

    ++maxTestLabel;

    ImageIteratorType refItr(refSegmentation, refSegmentation->GetLargestPossibleRegion());
    std::vector <unsigned int> labelsOverlap(maxTestLabel,0);
    std::vector <unsigned int> labelsNonOverlapping(maxTestLabel,0);
    std::vector <unsigned int> labelsSizes(maxTestLabel,0);

    ImageType::Pointer subImage = ImageType::New();
    subImage->Initialize();
    subImage->SetRegions(testSegmentation->GetLargestPossibleRegion());
    subImage->SetSpacing (testSegmentation->GetSpacing());
    subImage->SetOrigin (testSegmentation->GetOrigin());
    subImage->SetDirection (testSegmentation->GetDirection());
    subImage->Allocate();
    ImageIteratorType subItr(subImage,testSegmentation->GetLargestPossibleRegion());

    testItr.GoToBegin();
    while (!refItr.IsAtEnd())
    {
        if (refItr.Get() != 0)
            labelsOverlap[testItr.Get()]++;
        else
            labelsNonOverlapping[testItr.Get()]++;

        labelsSizes[testItr.Get()]++;

        ++refItr;
        ++testItr;
    }

    std::cout << "Processing " << maxTestLabel << " lesions in second timepoint..." << std::endl;
    for (unsigned int i = 1;i < maxTestLabel;++i)
    {
        testItr.GoToBegin();
        double ratioNonOverlapOverlap = static_cast <double> (labelsNonOverlapping[i]) / labelsOverlap[i];
        // Shrinking or not enough change -> label = maxTestLabel + 1
        if ((labelsNonOverlapping[i] <= minSizeInVoxel) || (ratioNonOverlapOverlap <= betaArg.getValue()))
        {
            while (!testItr.IsAtEnd())
            {
                if (testItr.Get() == i)
                    testItr.Set(maxTestLabel + 1);

                ++testItr;
            }

            continue;
        }

        // New lesion -> put it all as new -> label = maxTestLabel + 2
        double ratioOverlapSizes = static_cast <double> (labelsOverlap[i]) / labelsSizes[i];
        if (ratioOverlapSizes <= alphaArg.getValue())
        {
            while (!testItr.IsAtEnd())
            {
                if (testItr.Get() == i)
                    testItr.Set(maxTestLabel + 2);

                ++testItr;
            }

            continue;
        }

        // Test for growing lesion too large
        if (ratioNonOverlapOverlap > gammaArg.getValue())
        {
            // New lesion -> label = maxTestLabel + 2
            while (!testItr.IsAtEnd())
            {
                if (testItr.Get() == i)
                    testItr.Set(maxTestLabel + 2);

                ++testItr;
            }

            continue;
        }

        // We are looking at a growing lesion -> label = maxTestLabel + 3
        subImage->FillBuffer(0);
        refItr.GoToBegin();
        subItr.GoToBegin();
        while (!testItr.IsAtEnd())
        {
            if ((testItr.Get() == i)&&(refItr.Get() == 0))
                subItr.Set(1);

            ++testItr;
            ++refItr;
            ++subItr;
        }

        CCFilterType::Pointer tmpCCFilter = CCFilterType::New();
        tmpCCFilter->SetInput(subImage);
        tmpCCFilter->SetFullyConnected(fullConnectArg.isSet());
        tmpCCFilter->SetNumberOfThreads(nbpArg.getValue());
        tmpCCFilter->Update();

        // Remove too small test objects
        RelabelComponentFilterType::Pointer relabelTmpCCFilter = RelabelComponentFilterType::New();
        relabelTmpCCFilter->SetInput(tmpCCFilter->GetOutput());
        relabelTmpCCFilter->SetMinimumObjectSize(0);
        relabelTmpCCFilter->SetNumberOfThreads(nbpArg.getValue());
        relabelTmpCCFilter->Update();

        ImageIteratorType subCCItr(relabelTmpCCFilter->GetOutput(),testSegmentation->GetLargestPossibleRegion());
        std::vector <unsigned int> numVoxelsCCSub(1);
        while (!subCCItr.IsAtEnd())
        {
            unsigned int currentSize = numVoxelsCCSub.size();
            unsigned int value = subCCItr.Get();

            if (value >= currentSize)
            {
                numVoxelsCCSub.resize(value + 1);
                for (unsigned int j = currentSize;j <= value;++j)
                    numVoxelsCCSub[j] = 0;
            }

            numVoxelsCCSub[value]++;
            ++subCCItr;
        }

        unsigned int numSubLabels = numVoxelsCCSub.size() - 1;
        double ratioOverlap = static_cast <double> (numVoxelsCCSub[numSubLabels]) / labelsOverlap[i];
        bool okGrowing = (ratioOverlap > betaArg.getValue());

        refItr.GoToBegin();
        testItr.GoToBegin();
        if (okGrowing)
        {
            while (!refItr.IsAtEnd())
            {
                if (testItr.Get() == i)
                {
                    if (refItr.Get() == 0)
                        testItr.Set(maxTestLabel + 3);
                    else
                        testItr.Set(maxTestLabel + 1);
                }

                ++refItr;
                ++testItr;
            }
        }
        else
        {
            // Non growing lesion, setting to
            while (!testItr.IsAtEnd())
            {
                if (testItr.Get() == i)
                    testItr.Set(maxTestLabel + 1);

                ++testItr;
            }
        }
    }

    testItr.GoToBegin();
    while (!testItr.IsAtEnd())
    {
        unsigned int value = testItr.Get();
        if (value > 0)
            testItr.Set(value - maxTestLabel);

        ++testItr;
    }

    std::cout << "Writing output to " << outArg.getValue() << std::endl;

    anima::writeImage <ImageType> (outArg.getValue(),testSegmentation);
    return EXIT_SUCCESS;
}

