#include <tclap/CmdLine.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <iostream>
#include <fstream>

int main(int argc, char **argv)
{
  TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);
    
  TCLAP::ValueArg<std::string> inputArg("i","inputimage","Input image",true,"","Input image",cmd);
  TCLAP::ValueArg<std::string> outputArg("o","outputimage","Output image",true,"","Output image",cmd);
  TCLAP::ValueArg<std::string> maskArg("m","maskfile","mask file",false,"","mask file",cmd);

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
  typedef itk::ImageFileReader <DoubleImageType> DoubleReaderType;
  typedef itk::ImageFileReader <UCImageType> UCReaderType;
  typedef itk::ImageFileWriter <UCImageType> UCWriterType;
        
  typedef itk::OtsuThresholdImageFilter <DoubleImageType, UCImageType, UCImageType> OtsuThresholdImageFilterType;
	
  DoubleReaderType::Pointer inputImageReader = DoubleReaderType::New();
  inputImageReader->SetFileName(inputArg.getValue());

  try
    {
      inputImageReader->Update();
    }
  catch (itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      return 1;
    }

  UCReaderType::Pointer maskImageReader = UCReaderType::New();
  maskImageReader->SetFileName(maskArg.getValue());

  try
    {
      maskImageReader->Update();
    }
  catch (itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      return 1;
    }
  
  DoubleImageType::RegionType tmpRegionInputImage = inputImageReader->GetOutput()->GetLargestPossibleRegion();
  UCImageType::RegionType tmpRegionMaskImage = maskImageReader->GetOutput()->GetLargestPossibleRegion();

  if(tmpRegionInputImage.GetSize()[0]!=tmpRegionMaskImage.GetSize()[0] || tmpRegionInputImage.GetSize()[1]!=tmpRegionMaskImage.GetSize()[1] || tmpRegionInputImage.GetSize()[2]!=tmpRegionMaskImage.GetSize()[2])
    {
      std::cerr << "InputImage size != MaskImage size" << std::endl;
      return -1;
    }
  
  OtsuThresholdImageFilterType::Pointer otsuThrFilter = OtsuThresholdImageFilterType::New();
  otsuThrFilter->SetInput(inputImageReader->GetOutput());
  otsuThrFilter->SetMaskImage(maskImageReader->GetOutput());
	
  try
    {
      otsuThrFilter->Update();
    }
  catch (itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      return 1;
    }
      
  if (invArg.isSet())
    {
      UCImageIterator resIt(otsuThrFilter->GetOutput(),tmpRegionInputImage);
      unsigned int totalSize = tmpRegionInputImage.GetSize()[0]*tmpRegionInputImage.GetSize()[1]*tmpRegionInputImage.GetSize()[2];
      
      for (unsigned int i = 0;i < totalSize;++i)
        {
	  resIt.Set(1-resIt.Get());
	  ++resIt;
        }
    }
	
  UCWriterType::Pointer tmpWriter = UCWriterType::New();
	
  tmpWriter->SetInput(otsuThrFilter->GetOutput());
  tmpWriter->SetUseCompression(true);
  tmpWriter->SetFileName(outputArg.getValue());

  try
    {
      tmpWriter->Update();
    }
  catch (itk::ExceptionObject &e)
    {
      std::cerr << e << std::endl;
      return 1;
    }
  
  return 0;
}
