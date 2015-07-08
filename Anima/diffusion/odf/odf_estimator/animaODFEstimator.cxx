#include <animaODFEstimatorImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbe.h>
#include <tclap/CmdLine.h>

#include <animaGradientFileReader.h>
#include <animaReadWriteFunctions.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',"1.0");
	
    TCLAP::ValueArg<std::string> inArg("i","inputlist","File containing a list of diffusion weighted images",true,"","input diffusion images",cmd);
	TCLAP::ValueArg<std::string> gradArg("g","gradientlist","List of gradients (text file)",true,"","list of gradients",cmd);
	TCLAP::ValueArg<std::string> resArg("o","outputfile","Result ODF image",true,"","result ODF image",cmd);
    
    TCLAP::ValueArg<float> lambdaArg("l","lambda","Lambda regularization parameter (see Descoteaux MRM 2007)",false,0.006,"lambda for regularization",cmd);
    TCLAP::ValueArg<unsigned int> orderArg("k","order","Order of spherical harmonics basis",false,4,"Order of SH basis",cmd);
    
	TCLAP::ValueArg<double> sharpFactorArg("s","sharpenratio","Ratio for sharpening ODFs (see Descoteaux TMI 2009, default : 0.255)",false,0.255,"sharpening ratio",cmd);
	TCLAP::SwitchArg sharpenArg("S","sharpenodf","Sharpen ODF ? (default: no)",cmd,false);
	
	TCLAP::ValueArg<std::string> normSphereArg("n","normalizefile","Sphere tesselation for normalization",false,"","Normalization sphere file",cmd);
	TCLAP::SwitchArg normalizeArg("N","normalize","Normalize ODF ? (default: no)",cmd,false);
    
	TCLAP::SwitchArg radialArg("R","radialestimation","Use radial estimation (see Aganj et al) ? (default: no)",cmd,false);
	TCLAP::ValueArg<double> aganjRegFactorArg("d","adr","Delta threshold for signal regularization, only use if R option activated (see Aganj et al, default : 0.001)",false,0.001,"delta signal regularization",cmd);
	
    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default: all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
    
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }
    
    typedef anima::ODFEstimatorImageFilter <float,float> MainFilterType;
    typedef MainFilterType::TInputImage InputImageType;
	typedef itk::ImageFileWriter < MainFilterType::TOutputImage > ImageWriterType;
	
	MainFilterType::Pointer mainFilter = MainFilterType::New();
	mainFilter->SetLambda(lambdaArg.getValue());
	if (orderArg.getValue() % 2 == 0)
		mainFilter->SetLOrder(orderArg.getValue());
	else
		mainFilter->SetLOrder(orderArg.getValue() - 1);
	
	mainFilter->SetSharpen(sharpenArg.isSet());
	if (sharpenArg.isSet())
		mainFilter->SetSharpnessRatio(sharpFactorArg.getValue());
    
	mainFilter->SetUseAganjEstimation(radialArg.isSet());
    mainFilter->SetDeltaAganjRegularization(aganjRegFactorArg.getValue());
	
    mainFilter->SetNormalize(normalizeArg.isSet());
    if (normalizeArg.isSet())
        mainFilter->SetFileNameSphereTesselation(normSphereArg.getValue());
    
    typedef itk::Image<float,4> Image4DType;

    anima::setMultipleImageFilterInputsFromFileName<InputImageType,Image4DType,MainFilterType>(inArg.getValue(), mainFilter);
    
    typedef anima::GradientFileReader < std::vector < double >, double > GFReaderType;
    GFReaderType gfReader;
    gfReader.SetGradientFileName(gradArg.getValue());
    gfReader.SetGradientIndependentNormalization(true);
    
    gfReader.Update();
    
    GFReaderType::GradientVectorType directions = gfReader.GetGradients();
    
    for(unsigned int i = 0;i < directions.size();++i)
        mainFilter->AddGradientDirection(i,directions[i]);
    
	itk::TimeProbe tmpTime;
	tmpTime.Start();
	
	mainFilter->SetNumberOfThreads(nbpArg.getValue());
	mainFilter->Update();
	
	tmpTime.Stop();
	
	std::cout << "Execution Time: " << tmpTime.GetTotal() << std::endl;
	
	ImageWriterType::Pointer imageWriter = ImageWriterType::New();
	imageWriter->SetInput(mainFilter->GetOutput());
	imageWriter->SetFileName(resArg.getValue());
	imageWriter->SetUseCompression(true);
	
	imageWriter->Update();
	
	return 0;
}
