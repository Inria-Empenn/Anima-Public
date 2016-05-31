#include <iostream>
#include <tclap/CmdLine.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbe.h>
#include <itkImageConstIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include "animaReadWriteFunctions.h"

using namespace std;

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("Given a mask creates and a list of volume data, populates a csv file where each colmuns represents one of the volume data and each line a given voxel. The voxel for which the mask is 0 are not set into the csv file. \\ Note: the voxel location are lost into the csv file.", ' ',"1.0");
    TCLAP::ValueArg<std::string> dataListArg("i","database","Image List",true,"","image list",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","maskname","Mask",true,"","mask",cmd);
    TCLAP::ValueArg<std::string> resNameArg("o","outputname","CSV file name",true,"","file name for the output csv file",cmd);
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    std::string dataLTName = dataListArg.getValue();
    std::string maskName = maskArg.getValue();
    

    // load mask file
    typedef itk::Image <unsigned char, 3> MaskType;
    typedef itk::ImageFileReader < MaskType > itkMaskReader;
    itkMaskReader::Pointer maskRead = itkMaskReader::New();
    maskRead->SetFileName(maskName.c_str());
    maskRead->Update();
    typedef itk::ImageRegionConstIteratorWithIndex<MaskType>  MaskIterator;

    MaskIterator iterMaskImage( maskRead->GetOutput(), maskRead->GetOutput()->GetRequestedRegion() );
    

    typedef itk::Image <double, 3> ImageType;
    typedef itk::ImageRegionConstIterator<ImageType>  ImageIterator;

    // load filenames file
    ifstream fileIn(dataLTName.c_str());


    std::vector< std::vector<double> > allSamples;
    std::string oneLine;

    unsigned int nbImg=0;
    std::vector<std::string> labels;
    std::vector<unsigned int> xSample;  std::vector<unsigned int> ySample;  std::vector<unsigned int> zSample;
    while (std::getline(fileIn, oneLine))
    {
        std::cout<<nbImg<<std::endl;
        std::size_t pos = oneLine.find(":");               // position of ":" in str
        std::string imgeName = oneLine.substr (pos+1);     // get from ":" to the end
        labels.push_back(oneLine.substr(0,pos));


        itk::SmartPointer<ImageType> ltReader;
        try
        {
            ltReader=anima::readImage<ImageType>(imgeName);
        }
        catch(exception& e)
        {
            std::cout<<"file for field "<<labels.back()<<" is missing"<<std::endl;
            labels.pop_back();
            continue;
        }
        ImageIterator iterCurrentImage( ltReader, ltReader->GetRequestedRegion() );

        iterMaskImage.GoToBegin();
        iterCurrentImage.GoToBegin();
        std::vector<double> oneSample;
        while(!iterCurrentImage.IsAtEnd())
        {
            if(iterMaskImage.Value()!=0)
            {    oneSample.push_back(iterCurrentImage.Value());

                if(nbImg==0)
                {
                    const unsigned int x=iterCurrentImage.GetIndex()[0];
                    const unsigned int y=iterCurrentImage.GetIndex()[1];
                    const unsigned int z=iterCurrentImage.GetIndex()[2];
                    xSample.push_back(x);
                    ySample.push_back(y);
                    zSample.push_back(z);
                }
            }

            ++iterMaskImage;
            ++iterCurrentImage;
        }
        allSamples.push_back(oneSample);
        ++nbImg;
    }

    if (resNameArg.getValue() != "")
    {
        // save csv
        std::ofstream file(resNameArg.getValue().c_str());

        for(unsigned int indexI=0;indexI<allSamples.size();++indexI)
        {
            file<<labels[indexI];
            if(indexI<(allSamples.size()-1))
                file<<",";
        }
        file<<",indexX,indexY,indexZ";
        file<<std::endl;

        for(unsigned int indexJ=0;indexJ<allSamples[0].size();++indexJ)
        {
            for(unsigned int indexI=0;indexI<allSamples.size();++indexI)
            {

                file<<allSamples[indexI][indexJ];
                if(indexI<(allSamples.size()-1))
                    file<<",";
            }
            file<<","<<xSample[indexJ]<<","<<ySample[indexJ]<<","<<zSample[indexJ];
            file<<std::endl;
        }

        file.close();
    }

    return EXIT_SUCCESS;
}

