#include <animaReadWriteFunctions.h>

#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkMultiThreader.h>
#include <vnl/vnl_matrix.h>

#include <tclap/CmdLine.h>
#include <fstream>
#include <algorithm>

struct
{
    bool operator()(const std::pair <unsigned int, unsigned int> &a, const std::pair <unsigned int, unsigned int> &b) const
    {
        return a.second > b.second;
    }
} pair_decreasing_comparator;

int main(int argc, char * *argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ', ANIMA_VERSION);

    TCLAP::ValueArg <std::string> refArg("r", "ref", "Reference binary segmentation", true, "", "reference segmentation", cmd);
    TCLAP::ValueArg <std::string> testArg("t", "test", "Test binary segmentation", true, "", "test segmentation", cmd);
    TCLAP::ValueArg <std::string> outArg("o", "out-csv", "output CSV data on lesions", true, "", "output CSV", cmd);

    TCLAP::ValueArg <double> alphaArg("a", "alpha", "Alpha threshold (in between 0 and 1, default: 0.1)", false, 0.1, "alpha threshold", cmd);
    TCLAP::ValueArg <double> gammaArg("g", "gamma", "Gamma threshold (in between 0 and 1, default: 0.65)", false, 0.65, "gamma threshold", cmd);
    TCLAP::ValueArg <double> betaArg("b", "beta", "Beta threshold (in between 0 and 1, default: 0.7)", false, 0.7, "beta threshold", cmd);

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
    typedef itk::ImageRegionConstIterator <ImageType> ImageIteratorType;

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
    unsigned int minSizeInVoxel = (unsigned int)std::ceil(minVolumeArg.getValue() / spacingTot);

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

    ImageIteratorType refItr(refSegmentation, refSegmentation->GetLargestPossibleRegion());
    ImageIteratorType testItr(testSegmentation, testSegmentation->GetLargestPossibleRegion());

    unsigned int maxRefLabel = 0;
    unsigned int maxTestLabel = 0;
    while (!refItr.IsAtEnd())
    {
        if (refItr.Get() > maxRefLabel)
            maxRefLabel = refItr.Get();

        if (testItr.Get() > maxTestLabel)
            maxTestLabel = testItr.Get();

        ++refItr;
        ++testItr;
    }

    ++maxRefLabel;
    ++maxTestLabel;
    if (maxRefLabel <= 1)
        return EXIT_FAILURE;

    refItr.GoToBegin();
    testItr.GoToBegin();

    vnl_matrix <unsigned int> labelsOverlap(maxRefLabel,maxTestLabel);
    labelsOverlap.fill(0);

    while (!refItr.IsAtEnd())
    {
        labelsOverlap(refItr.Get(),testItr.Get())++;

        ++refItr;
        ++testItr;
    }

    std::vector < std::pair <unsigned int, unsigned int> > subVectorForSort(maxTestLabel - 1);
    std::vector < std::pair <double,unsigned int> > detectionTable(maxRefLabel-1);
    std::vector <double> columnSums(maxTestLabel,0);
    std::vector <double> lineSums(maxRefLabel,0);

    for (unsigned int i = 1;i < maxRefLabel;++i)
    {
        for (unsigned int j = 1;j < maxTestLabel;++j)
        {
            lineSums[i] += labelsOverlap(i,j);
            columnSums[j] += labelsOverlap(i,j);
        }
    }

    for (unsigned int i = 1;i < maxRefLabel;++i)
    {
        double denom = labelsOverlap(i,0) + lineSums[i];
        double Si = lineSums[i] / denom;

        if (Si <= alphaArg.getValue())
        {
            detectionTable[i-1] = std::make_pair(denom * spacingTot,0);
            continue;
        }

        for (unsigned int j = 1;j < maxTestLabel;++j)
            subVectorForSort[j-1] = std::pair <unsigned int,unsigned int>(j,labelsOverlap(i,j));

        std::sort(subVectorForSort.begin(),subVectorForSort.end(),pair_decreasing_comparator);

        double wSum = 0;
        unsigned int detectedObject = 1;
        unsigned int k = 0;
        while (wSum < gammaArg.getValue())
        {
            unsigned int kIndex = subVectorForSort[k].first;
            double sumK = labelsOverlap(0,kIndex) + columnSums[kIndex];
            double Tk = labelsOverlap(0,kIndex) / sumK;

            if (Tk > betaArg.getValue())
            {
                detectedObject = 0;
                break;
            }

            wSum += labelsOverlap(i,kIndex) / lineSums[i];
            ++k;
        }

        detectionTable[i-1] = std::make_pair(denom * spacingTot,detectedObject);
    }

    unsigned int totalNumberOfDetections = 0;
    unsigned int refCount = 0;
    for (unsigned int i = 1;i < maxRefLabel;++i)
    {
        if (detectionTable[i].second > 0)
            ++totalNumberOfDetections;

        for (unsigned int j = 0;j < maxTestLabel;++j)
            refCount += labelsOverlap(i,j);
    }

    std::ofstream outputFile(outArg.getValue());
    outputFile.precision(10);
    outputFile << "Total number of reference objects:," << maxRefLabel - 1 << std::endl;
    outputFile << "Total number of detections:," << totalNumberOfDetections << std::endl;
    outputFile << "Total objects volume:," << refCount * spacingTot << std::endl;

    outputFile << "Object number,Object volume,Detected flag" << std::endl;
    for (unsigned int i = 1;i < maxRefLabel;++i)
        outputFile << i << "," << detectionTable[i-1].first << "," << detectionTable[i-1].second << std::endl;

    outputFile.close();

    return EXIT_SUCCESS;
}

