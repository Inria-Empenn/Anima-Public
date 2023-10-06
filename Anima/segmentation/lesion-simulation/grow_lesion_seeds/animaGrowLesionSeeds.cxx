#include <animaInhomogeneousDiffusionImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inArg("i","inputfile","Input seed image",true,"","input seed image",cmd);
    TCLAP::ValueArg<std::string> probaArg("p","probafile","Input probability image",true,"","input probability image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","Output lesions image",true,"","output lesions image",cmd);
	
    TCLAP::ValueArg<unsigned short> numStepsArg("s","numsteps","Number of steps (default: 100)",false,100,"number of steps",cmd);
    TCLAP::ValueArg<double> stepLengthArg("l","steplength","Time step (default: 0.1)",false,0.1,"time step",cmd);
    TCLAP::ValueArg<double> diffSourceArg("d","diffsource","Source of diffusion intensity (default: 0)",false,0,"diffusion source intensity",cmd);
    
    TCLAP::ValueArg<unsigned int> nbpArg("T","nthreads","Number of cores to run on (default: all cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"Number of cores",cmd);
	
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }
	
	typedef itk::Image <unsigned short,3> ImageType;
    typedef itk::ImageFileReader <ImageType> ImageReaderType;
        
    ImageReaderType::Pointer inputRead = ImageReaderType::New();
    inputRead->SetFileName(inArg.getValue());
    inputRead->Update();
        
	typedef itk::Image <double,3> ProbaImageType;
    typedef itk::ImageFileReader <ProbaImageType> ProbaImageReaderType;
    typedef itk::ImageFileWriter <ProbaImageType> ProbaImageWriterType;
    
    ProbaImageReaderType::Pointer probaRead = ProbaImageReaderType::New();
    probaRead->SetFileName(probaArg.getValue());
    probaRead->Update();

    typedef anima::InhomogeneousDiffusionImageFilter <ImageType, ProbaImageType, ProbaImageType> MainFilterType;
    
    MainFilterType::Pointer mainFilter = MainFilterType::New();
    
    mainFilter->SetInput(inputRead->GetOutput());
    mainFilter->SetNumberOfWorkUnits(nbpArg.getValue());
    mainFilter->SetNumberOfSteps(numStepsArg.getValue());
    mainFilter->SetStepLength(stepLengthArg.getValue());
    mainFilter->SetDiffusionSourceFactor(diffSourceArg.getValue());
    mainFilter->SetDiffusionScalarsImage(probaRead->GetOutput());
    
    mainFilter->Update();
    
    ProbaImageWriterType::Pointer diffWriter = ProbaImageWriterType::New();
    diffWriter->SetInput(mainFilter->GetOutput());
    diffWriter->SetFileName(outArg.getValue());
    diffWriter->SetUseCompression(true);

    diffWriter->Update();
	
	return 0;
}
