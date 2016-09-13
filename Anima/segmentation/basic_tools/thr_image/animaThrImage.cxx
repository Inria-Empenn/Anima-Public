#include <itkThresholdImageFilter.h>
#include <itkThresholdLabelerImageFilter.h>
#include <itkImageRegionConstIterator.h>

#include <animaReadWriteFunctions.h>

#include <limits.h>
#include <iostream>
#include <fstream>
#include <tclap/CmdLine.h>

using namespace std;

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inputArg("i","inputimage","Input image",true,"","Input image",cmd);
    TCLAP::ValueArg<std::string> outputArg("o","outputimage","Output image",true,"","Output image",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","maskfile","mask file",false,"","mask file",cmd);

    TCLAP::ValueArg<double> thrArg("t","thr","Threshold value",false,1.0,"Threshold value",cmd);
    TCLAP::ValueArg<double> upperThrArg("u","uthr","Upper threshold value",false,USHRT_MAX,"Upper threshold value",cmd);
    TCLAP::ValueArg<double> adaptThrArg("a","adaptivethr","Adaptative threshold value (between 0 and 1)",false,0.0,"adaptative threshold value",cmd);
    
    TCLAP::SwitchArg invArg("I","inv","Computes 1-res",cmd,false);
    
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    typedef itk::Image<double,3> DoubleImageType;
    typedef itk::Image<unsigned char, 3> UCImageType;
    
    typedef itk::ImageRegionConstIterator <DoubleImageType> DoubleImageIterator;
    typedef itk::ImageRegionIterator <UCImageType> UCImageIterator;

    typedef itk::ThresholdImageFilter <DoubleImageType> ThresholdFilterType;
    typedef itk::ThresholdLabelerImageFilter <DoubleImageType,UCImageType> LabelerFilterType;

    DoubleImageType::Pointer inputImage = anima::readImage <DoubleImageType> (inArg.getValue());
    
    DoubleImageType::RegionType tmpRegionInputImage = inputImage->GetLargestPossibleRegion();
    DoubleImageIterator doubleIt(inputImage,tmpRegionInputImage);

    unsigned int totalSize = tmpRegionInputImage.GetSize()[0]*tmpRegionInputImage.GetSize()[1]*tmpRegionInputImage.GetSize()[2];
    std::vector <double> tmpVec (totalSize);

    unsigned int idx=0;
    
    if(maskArg.getValue()!="")
    {
	UCImageType::Pointer maskImage = anima::readImage <UCImageType> (maskArg.getValue());

        UCImageType::RegionType tmpRegionMaskImage = maskImage->GetLargestPossibleRegion();
        UCImageIterator ucIt(maskImage,tmpRegionMaskImage);

        if(tmpRegionInputImage.GetSize()[0]!=tmpRegionMaskImage.GetSize()[0] || tmpRegionInputImage.GetSize()[1]!=tmpRegionMaskImage.GetSize()[1] || tmpRegionInputImage.GetSize()[2]!=tmpRegionMaskImage.GetSize()[2])
        {
            std::cerr << "InputImage size != MaskImage size" << std::endl;
            return -1;
        }

        for (unsigned int i = 0;i < totalSize;++i)
        {
            if(ucIt.Get()==1)
            {
                tmpVec[idx] = doubleIt.Get();
                ++idx;
            }
            ++ucIt;
            ++doubleIt;
        }

        tmpVec.resize(idx);
    }
    else
    {
        idx=totalSize;
        for (unsigned int i = 0;i < totalSize;++i)
        {
            tmpVec[i] = doubleIt.Get();
            ++doubleIt;
        }
    }

    if(adaptThrArg.getValue() < 0 || adaptThrArg.getValue() > 1)
    {
        std::cerr << "Adaptative threshold value has to be included in the [0,1] interval" << std::endl;
        return -1;
    }
    
    unsigned int partialElt = (unsigned int)floor(adaptThrArg.getValue()*idx);
    if (partialElt == idx)
        partialElt = idx - 1;

    if (partialElt != 0)
        std::partial_sort(tmpVec.begin(),tmpVec.begin() + partialElt + 1,tmpVec.end());

    double thrV = thrArg.getValue();
    if (partialElt != 0)
        thrV = (thrV < tmpVec[partialElt]) ? tmpVec[partialElt] : thrV;

    ThresholdFilterType::Pointer thrFilter = ThresholdFilterType::New();
    thrFilter->SetInput(inputImage);

    double upperThr = upperThrArg.getValue();
    if (upperThr < thrV)
        upperThr = thrV;

    if (upperThr != USHRT_MAX)
        thrFilter->ThresholdOutside(thrV,upperThr);
    else
        thrFilter->ThresholdBelow(thrV);
    
    try
    {
        thrFilter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return 1;
    }
    
    LabelerFilterType::Pointer mainFilter = LabelerFilterType::New();
    mainFilter->SetInput(thrFilter->GetOutput());

    LabelerFilterType::RealThresholdVector thrVals;
    thrVals.push_back(0);

    mainFilter->SetRealThresholds(thrVals);

    try
    {
        mainFilter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return 1;
    }
    
    if (invArg.isSet())
    {
        UCImageIterator resIt(mainFilter->GetOutput(),tmpRegionInputImage);
        
        for (unsigned int i = 0;i < totalSize;++i)
        {
            resIt.Set(1-resIt.Get());
            ++resIt;
        }
    }

    anima::writeImage <UCImageType> (outArg.getValue(),mainFilter->GetOutput());
        
    return 0;
}
