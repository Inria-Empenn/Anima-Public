#include <animaQMRISampleCreationImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <fstream>

#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inArg("i","inputfile","Input QMRI images in text file",true,"","input QMRI images",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputprefix","Output QMR images prefix",true,"","output QMRI prefix",cmd);
    TCLAP::ValueArg<std::string> outLesionMaskArg("O","outputmask","Output lesion mask",false,"","output lesion mask",cmd);

    TCLAP::ValueArg<std::string> probaArg("p","probafile","Lesions probability image",true,"","lesions probability image",cmd);
    TCLAP::ValueArg<std::string> varImagesArg("v","varfiles","Input QMR variance images in a text file",true,"","input QMR variance images",cmd);

    TCLAP::ValueArg<std::string> lesionSizeDistArg("l","lesionsizedist","Text file with the cumulative distribution of lesion sizes",true,"","lesion size distribution",cmd);
    TCLAP::ValueArg<std::string> qmriLesionRelationshipsArg("r","lesionqmrirel","Linear relationships between QMRI intensities and distance to lesion border",true,"","qMRI / lesion relationships",cmd);
    
    TCLAP::ValueArg<unsigned int> numSeedsArg("n","numseeds","Number of lesion seeds (default: 5)",false,5,"Number of lesion seeds",cmd);
    TCLAP::ValueArg<unsigned int> lesionMinSizeArg("m","lesionminsize","Minimal number of voxels in a lesion (default: 5)",false,5,"minimal lesion size in voxels",cmd);
    TCLAP::ValueArg<double> minDistArg("d","dist","Minimal distance between seeds (default: 10 mm)",false,10,"minimal distance between seeds",cmd);

    TCLAP::ValueArg<double> diffThresholdArg("t","diffthr","Diffusion grower threshold to binarize grown lesion (default: 1)",false,1,"diffusion grower threshold",cmd);

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
	
	typedef itk::Image <double,3> ImageType;
    typedef itk::ImageFileReader <ImageType> ImageReaderType;
    
    typedef anima::QMRISampleCreationImageFilter < itk::Image <double, 3>,itk::Image <double, 3> > MainFilterType;
    MainFilterType::Pointer mainFilter = MainFilterType::New();
    
    unsigned int pos = 0;
    std::ifstream inputFile(inArg.getValue().c_str());
    while (!inputFile.eof())
    {
        char tmpStr[2048];
        inputFile.getline(tmpStr,2048);
        
        if (strcmp(tmpStr,"") != 0)
        {
            ImageReaderType::Pointer tmpReader = ImageReaderType::New();
            tmpReader->SetFileName(tmpStr);
            tmpReader->Update();
            
            std::cout << "Loading " << pos+1 << "th qMRI average image " << tmpStr << std::endl;
            mainFilter->SetInput(pos,tmpReader->GetOutput());
            ++pos;
        }
    }
    
    inputFile.close();
    pos = 0;
    
    std::ifstream varFile(varImagesArg.getValue().c_str());
    while (!varFile.eof())
    {
        char tmpStr[2048];
        varFile.getline(tmpStr,2048);
        
        if (strcmp(tmpStr,"") != 0)
        {
            ImageReaderType::Pointer tmpReader = ImageReaderType::New();
            tmpReader->SetFileName(tmpStr);
            tmpReader->Update();
            
            std::cout << "Loading " << pos+1 << "th qMRI variance image " << tmpStr << std::endl;
            mainFilter->AddQMRIVarianceImage(tmpReader->GetOutput());
            ++pos;
        }
    }
    
    varFile.close();
    
    mainFilter->SetQMRILesionRelationshipsFile(qmriLesionRelationshipsArg.getValue());
    mainFilter->ReadLesionSizesDistributions(lesionSizeDistArg.getValue());
    
    ImageReaderType::Pointer probaRead = ImageReaderType::New();
    probaRead->SetFileName(probaArg.getValue());
    probaRead->Update();
    
    mainFilter->SetLesionsProbabilityMap(probaRead->GetOutput());
    mainFilter->SetNumberOfWorkUnits(nbpArg.getValue());
    mainFilter->SetLesionDiffusionThreshold(diffThresholdArg.getValue());
    mainFilter->SetLesionMinimalSize(lesionMinSizeArg.getValue());
    mainFilter->SetMinimalDistanceBetweenLesions(minDistArg.getValue());
    mainFilter->SetNumberOfSeeds(numSeedsArg.getValue());
    
    mainFilter->Update();
    
    typedef itk::ImageFileWriter <ImageType> ImageWriterType;
    std::string outFileBaseName = outArg.getValue() + "_";
    
    for (unsigned int i = 0;i < mainFilter->GetNumberOfOutputs();++i)
    {
        std::ostringstream outNum;
        outNum << i;
        std::string outFileName = outFileBaseName;
        outFileName += outNum.str();
        outFileName += ".nii.gz";
        
        ImageWriterType::Pointer tmpWriter = ImageWriterType::New();
        tmpWriter->SetFileName(outFileName);
        tmpWriter->SetInput(mainFilter->GetOutput(i));
        tmpWriter->Update();
    }

	if (outLesionMaskArg.getValue() != "")
    {
        typedef itk::ImageFileWriter <MainFilterType::MaskImageType> MaskImageWriterType;
        
        MaskImageWriterType::Pointer maskWriter = MaskImageWriterType::New();
        maskWriter->SetInput(mainFilter->GetLesionsOutputMask());
        maskWriter->SetFileName(outLesionMaskArg.getValue());
        
        maskWriter->Update();
    }
    
	return 0;
}
