#include <iostream>
#include <algorithm>
#include <tclap/CmdLine.h>

#include <itkImageRegionConstIterator.h>
#include <vnl/vnl_matrix.h>
#include <animaReadWriteFunctions.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("Given a mask creates and a list of volume data, populates a csv file where each column represents one of the volume data and each line a given ROI. The values reported are selected from mean, variance, median."
                       " INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    TCLAP::ValueArg <std::string> dataListArg("i","input","image or image list",true,"","images",cmd);
    TCLAP::ValueArg <std::string> roiArg("r","roi","ROI image",true,"","ROI image",cmd);
    TCLAP::ValueArg <std::string> statArg("s","stat","Type of statistic computed (choose from median, [average], variance, max)",false,"average","statistic type",cmd);
    TCLAP::ValueArg <std::string> resNameArg("o","output","CSV file name",true,"","file name for the output csv file",cmd);
    TCLAP::ValueArg <unsigned int> precisionArg("p","precision","Precision of values output (integer, default: 12)",false,12,"precision",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    // load mask file
    using ROIImageType = itk::Image <unsigned short, 3>;
    ROIImageType::Pointer roiImage = anima::readImage <ROIImageType> (roiArg.getValue());
    using ROIIterator = itk::ImageRegionConstIterator <ROIImageType>;

    ROIIterator iterROIImage(roiImage,roiImage->GetRequestedRegion());

    std::vector <unsigned short> roiLabels;

    // First find labels in ROI image
    while (!iterROIImage.IsAtEnd())
    {
        unsigned short roiValue = iterROIImage.Get();
        if (roiValue == 0)
        {
            ++iterROIImage;
            continue;
        }

        bool foundLabel = false;
        for (unsigned int i = 0;i < roiLabels.size();++i)
        {
            if (roiValue == roiLabels[i])
            {
                foundLabel = true;
                break;
            }
        }

        if (!foundLabel)
            roiLabels.push_back(roiValue);

        ++iterROIImage;
    }

    std::sort(roiLabels.begin(),roiLabels.end());
    unsigned int numLabels = roiLabels.size();
    std::map <unsigned short, unsigned short> reverseLabels;
    for (unsigned int i = 0;i < numLabels;++i)
        reverseLabels[roiLabels[i]] = i;

    using ImageType = itk::Image <double, 3>;
    using ImageIterator = itk::ImageRegionConstIterator<ImageType>;

    // load filenames file    
    std::vector <std::string> dataFileNames;
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(dataListArg.getValue().c_str(), itk::IOFileModeEnum::ReadMode);

    if (imageIO)
        dataFileNames.push_back(dataListArg.getValue());
    else
    {
        std::ifstream fileIn(dataListArg.getValue());
        std::string oneLine;

        while (std::getline(fileIn, oneLine))
            dataFileNames.push_back(oneLine);
    }

    //Now process voxels for each image
    unsigned int numDataFiles = dataFileNames.size();
    std::vector < std::vector <double> > roiLabelValues(numLabels);
    vnl_matrix <double> statsValues(numLabels, numDataFiles);

    for (unsigned int i = 0;i < numDataFiles;++i)
    {
        for (unsigned int j = 0;j < numLabels;++j)
            roiLabelValues[j].clear();

        ImageType::Pointer currentImage = anima::readImage <ImageType> (dataFileNames[i]);
        ImageIterator imageIt(currentImage, currentImage->GetLargestPossibleRegion());

        iterROIImage.GoToBegin();
        while (!iterROIImage.IsAtEnd())
        {
            unsigned short roiValue = iterROIImage.Get();

            if (roiValue == 0)
            {
                ++iterROIImage;
                ++imageIt;

                continue;
            }

            roiLabelValues[reverseLabels[roiValue]].push_back(imageIt.Get());

            ++iterROIImage;
            ++imageIt;
        }

        // Now we have the data, let's compute the stats
        for (unsigned int j = 0;j < numLabels;++j)
        {
            if (statArg.getValue() == "median")
            {
                unsigned int numValues = roiLabelValues[j].size();
                unsigned int indexTaken = std::floor(numValues / 2.0);
                std::partial_sort(roiLabelValues[j].begin(),roiLabelValues[j].begin() + indexTaken + 1,roiLabelValues[j].end());
                statsValues(j,i) = roiLabelValues[j][indexTaken];
            }
            else if (statArg.getValue() == "max")
            {
                unsigned int numValues = roiLabelValues[j].size();
                double maxValue = roiLabelValues[j][0];
                for (unsigned int k = 1;k < numValues;++k)
                {
                    if (maxValue < roiLabelValues[j][k])
                        maxValue = roiLabelValues[j][k];
                }
            }
            else
            {
                double meanValue = 0.0;
                unsigned int numValues = roiLabelValues[j].size();
                for (unsigned int k = 0;k < numValues;++k)
                    meanValue += roiLabelValues[j][k];

                meanValue /= numValues;

                if (statArg.getValue() == "variance")
                {
                    double varianceValue = 0.0;
                    for (unsigned int k = 0;k < numValues;++k)
                        varianceValue += (meanValue - roiLabelValues[j][k]) * (meanValue - roiLabelValues[j][k]);

                    if (numValues > 1)
                        varianceValue /= (numValues - 1.0);

                    statsValues(j,i) = varianceValue;
                }
                else
                    statsValues(j,i) = meanValue;
            }
        }
    }

    // save csv
    std::ofstream file(resNameArg.getValue());
    file.precision(precisionArg.getValue());

    file << "Label,";
    for (unsigned int i = 0;i < numDataFiles;++i)
    {
        file << dataFileNames[i];
        if (i != numDataFiles - 1)
            file << ",";
        else
            file << std::endl;
    }

    for (unsigned int i = 0;i < numLabels;++i)
    {
        file << roiLabels[i] << ",";
        for (unsigned int j = 0;j < numDataFiles;++j)
        {
            file << statsValues(i,j);
            if (j != numDataFiles - 1)
                file << ",";
            else
                file << std::endl;
        }
    }

    file.close();

    return EXIT_SUCCESS;
}

