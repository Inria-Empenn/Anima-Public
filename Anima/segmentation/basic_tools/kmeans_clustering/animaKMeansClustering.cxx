#include <tclap/CmdLine.h>

#include <itkScalarImageKmeansImageFilter.h>
#include <animaReadWriteFunctions.h>
#include <itkImageRegionIterator.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inArg("i","inputfile","Input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","output image",true,"","output image",cmd);
    TCLAP::ValueArg<unsigned int> numClassArg("c","numclass","Number of classes",false,4,"number of classes",cmd);
    TCLAP::ValueArg<int> keptClassArg("k","class","Class to keep",false,-1,"class to be kept",cmd);
    TCLAP::ValueArg<unsigned int> nbpArg("n","ncores","Number of cores (default: all cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of cores",cmd);
    
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }
    
    std::string refName, resName;
    refName = inArg.getValue();
    resName = outArg.getValue();
    
    typedef itk::Image <float,3> FloatImageType;
    typedef itk::ScalarImageKmeansImageFilter <FloatImageType> MainFilterType;
    
    MainFilterType::Pointer kmeansFilter = MainFilterType::New();
    kmeansFilter->SetInput(anima::readImage <FloatImageType> (refName));
    kmeansFilter->SetUseNonContiguousLabels(false);

    for (unsigned int i = 0;i < numClassArg.getValue();++i)
        kmeansFilter->AddClassWithInitialMean(i * 100);

    kmeansFilter->SetNumberOfWorkUnits(nbpArg.getValue());
    kmeansFilter->Update();
    
    MainFilterType::ParametersType estimatedMeans = kmeansFilter->GetFinalMeans();
    typedef MainFilterType::OutputImageType OutputImageType;

    if (keptClassArg.getValue() >= 0)
    {
        const unsigned int numberOfClasses = estimatedMeans.Size();
        unsigned int numLabel = 0;
    
        for(unsigned int i = 0 ; i < numberOfClasses ; ++i)
        {
            unsigned int numInf = 0;
        
            for (unsigned int j = 0;j < numberOfClasses;++j)
            {
                if (estimatedMeans[j] < estimatedMeans[i])
                    ++numInf;
            }
        
            if (numInf == keptClassArg.getValue())
            {
                numLabel = i;
                break;
            }
        }

        typedef itk::ImageRegionIterator <OutputImageType> IteratorType;
        IteratorType outIterator(kmeansFilter->GetOutput(),kmeansFilter->GetOutput()->GetLargestPossibleRegion());
    
        while(!outIterator.IsAtEnd())
        {
            outIterator.Set((outIterator.Get() == numLabel));
            ++outIterator;
        }
    }

    anima::writeImage <OutputImageType> (resName,kmeansFilter->GetOutput());

    return EXIT_SUCCESS;
}
