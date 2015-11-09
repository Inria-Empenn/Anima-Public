#include <itkThresholdImageFilter.h>
#include <itkThresholdLabelerImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIterator.h>

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
	
	typedef itk::Image<double,3> InputImageType;
	typedef itk::Image<unsigned char, 3> OutputImageType;
	typedef itk::ImageFileReader <InputImageType> InputReaderType;
	typedef itk::ImageFileWriter <OutputImageType> OutputWriterType;
	typedef itk::ImageRegionConstIterator <InputImageType> InputImageIterator;
	typedef itk::ImageRegionIterator <OutputImageType> OutputImageIterator;
	
	typedef itk::ThresholdImageFilter <InputImageType> ThresholdFilterType;
	typedef itk::ThresholdLabelerImageFilter <InputImageType,OutputImageType> LabelerFilterType;
	
	
	InputReaderType::Pointer tmpRead = InputReaderType::New();
	tmpRead->SetFileName(inputArg.getValue());

    try
    {
        tmpRead->Update();
	}
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return 1;
    }
	
	InputImageType::RegionType tmpRegion = tmpRead->GetOutput()->GetLargestPossibleRegion();
	unsigned int totalSize = tmpRegion.GetSize()[0]*tmpRegion.GetSize()[1]*tmpRegion.GetSize()[2];
	std::vector <double> tmpVec (totalSize);
	
	InputImageIterator tmpIt(tmpRead->GetOutput(),tmpRegion);
	for (unsigned int i = 0;i < totalSize;++i)
	{
		tmpVec[i] = tmpIt.Get();
		++tmpIt;
	}
	
	unsigned int partialElt = (unsigned int)floor(adaptThrArg.getValue()*totalSize);
	if (partialElt == totalSize)
		partialElt = totalSize - 1;
	
	if (partialElt != 0)
		std::partial_sort(tmpVec.begin(),tmpVec.begin() + partialElt + 1,tmpVec.end());
	
	double thrV = thrArg.getValue();
	if (partialElt != 0)
		thrV = (thrV < tmpVec[partialElt]) ? tmpVec[partialElt] : thrV;
	
	ThresholdFilterType::Pointer thrFilter = ThresholdFilterType::New();
	thrFilter->SetInput(tmpRead->GetOutput());
	
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
        OutputImageIterator resIt(mainFilter->GetOutput(),tmpRegion);
        
        for (unsigned int i = 0;i < totalSize;++i)
        {
            resIt.Set(1-resIt.Get());
            ++resIt;
        }
    }
	
	OutputWriterType::Pointer tmpWriter = OutputWriterType::New();
	
	tmpWriter->SetInput(mainFilter->GetOutput());
	tmpWriter->SetUseCompression(true);
	tmpWriter->SetFileName(outputArg.getValue());
	
	tmpWriter->Update();
    
	return 0;
}
