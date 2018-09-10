#include <cmath>

#include <tclap/CmdLine.h>

#include <itkVectorImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <animaStimulatedSpinEchoImageFilter.h>

int main(int argc, char *argv[] )
{
    TCLAP::CmdLine cmd("Simulator for stimulated spin echo T2 sequences\nINRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> t1MapArg("","t1","Input T1 map",true,"","T1 map",cmd);
    TCLAP::ValueArg<std::string> t2MapArg("","t2","Input T2 map",true,"","T2 map",cmd);
    TCLAP::ValueArg<std::string> m0ImageArg("","m0","Input M0 image",true,"","M0 image",cmd);
    TCLAP::ValueArg<std::string> b1ImageArg("","b1","Input B1 image",false,"","B1 image",cmd);
    TCLAP::ValueArg<std::string> resArg("o","output","Output simulated image",true,"","output simulated image",cmd);
    
    TCLAP::ValueArg<unsigned int> neArg("n","ne","Number of echoes, default: 7",false, 7,"number of echoes",cmd);
    TCLAP::ValueArg<double> teArg("e","te","TE echo spacing (ms), default: 10ms",false, 10,"TE value",cmd);
    TCLAP::ValueArg<double> xfaArg("x","excite-fa","Excitation flip angle (degrees), default: 90 degrees", false, 90,"excitation flip angle value",cmd);
    TCLAP::ValueArg<double> faArg("f","fa","Flip angle (degrees), default: 180 degrees", false, 180,"flip angle value",cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("T","numberofthreads","Number of threads to run on (default : all cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }
    
    // Parse arguments
    std::string input1FileName = t1MapArg.getValue();
    std::string input2FileName = t2MapArg.getValue();
    std::string input3FileName = m0ImageArg.getValue();
    std::string outputFileName = resArg.getValue();
    
    // Setup types
    typedef itk::Image<double, 3> ImageType;
    typedef itk::VectorImage<double, 3> OutputImageType;
    typedef anima::StimulatedSpinEchoImageFilter<ImageType,OutputImageType> FilterType;
    
    // Read T1 map
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer reader1 = ReaderType::New();
    reader1->SetFileName(input1FileName);
    reader1->Update();
    
    // Read T2s map
    ReaderType::Pointer reader2 = ReaderType::New();
    reader2->SetFileName(input2FileName);
    reader2->Update();
    
    // Read M0 map
    ReaderType::Pointer reader3 = ReaderType::New();
    reader3->SetFileName(input3FileName);
    reader3->Update();
    
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(0,reader1->GetOutput());
    filter->SetInput(1,reader2->GetOutput());
    filter->SetInput(2,reader3->GetOutput());
    filter->SetNumberOfWorkUnits(nbpArg.getValue());
    
    // Read B1 map
    if (b1ImageArg.getValue() != "")
    {
        ReaderType::Pointer readerB1 = ReaderType::New();
        readerB1->SetFileName(b1ImageArg.getValue().c_str());
        readerB1->Update();
        
        filter->SetInput(3,readerB1->GetOutput());
    }

    filter->SetEchoSpacing(teArg.getValue());
    filter->SetFlipAngle(faArg.getValue() * M_PI / 180.0);
    filter->SetNumberOfEchoes(neArg.getValue());
    filter->SetExcitationFlipAngle(xfaArg.getValue() * M_PI / 180.0);
    
    filter->Update();
    
    // Save the output image
    typedef itk::ImageFileWriter <FilterType::Image4DType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outputFileName);
    writer->SetInput(filter->GetOutputAs4DImage());
    writer->Update();
    
    // return
    return EXIT_SUCCESS;
}

