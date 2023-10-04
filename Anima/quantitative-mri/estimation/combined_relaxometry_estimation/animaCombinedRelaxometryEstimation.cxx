#include <cmath>

#include <iostream>
#include <tclap/CmdLine.h>

#include <animaCombinedRelaxometryEstimationImageFilter.h>
#include <animaReadWriteFunctions.h>
#include <itkTimeProbe.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> t1Arg("","t1","List of T1 relaxometry images",true,"","T1 relaxometry images",cmd);
    TCLAP::ValueArg<std::string> t2Arg("","t2","List of T2 relaxometry images",true,"","T2 relaxometry images",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask",false,"","computation mask",cmd);

    TCLAP::ValueArg<std::string> resT1Arg("","out-t1","Result T1 image",true,"","result T1 image",cmd);
    TCLAP::ValueArg<std::string> resT2Arg("","out-t2","Result T2 image",false,"","result T2 image",cmd);
    TCLAP::ValueArg<std::string> resT1M0Arg("","out-m0-t1","Result M0 image (for T1 image)",false,"","result T1 M0 image",cmd);
    TCLAP::ValueArg<std::string> resT2M0Arg("","out-m0-t2","Result M0 image (for T2 image)",false,"","result T2 M0 image",cmd);
    TCLAP::ValueArg<std::string> resB1Arg("","out-b1","Result B1 image",false,"","result B1 image",cmd);
    TCLAP::ValueArg<std::string> resB1AddArg("","out-b1-add","Result B1 additive image",false,"","result B1 additive image",cmd);

    TCLAP::ValueArg<double> trT1Arg("","tr-t1","Repetition time for T1 relaxometry (default: 5000)",false,5000,"T1 repetition time",cmd);
    TCLAP::ValueArg<std::string> t1FlipAnglesArg("","t1-flip","Text file listing flip angles for T1 (in degrees)",true,"","T1 flip angles list",cmd);

    TCLAP::ValueArg<double> trT2Arg("","tr-t2","Repetition time for T2 relaxometry (default: 5000)",false,5000,"T2 repetition time",cmd);
    TCLAP::ValueArg<double> t2EchoSpacingArg("e","t2-echo-spacing","Echo spacing between T2 relaxometry images (default: 50)",false,50,"T2 echo spacing",cmd);
    TCLAP::ValueArg<double> excitationT2FlipAngleArg("","t2-ex-flip","Excitation flip angle for T2 (in degrees, default: 90)",false,90,"T2 excitation flip angle",cmd);
    TCLAP::ValueArg<double> t2FlipAngleArg("","t2-flip","All flip angles for T2 (in degrees, default: 180)",false,180,"T2 flip angle",cmd);

    TCLAP::ValueArg<double> upperBoundT1Arg("","upper-bound-t1","T1 value upper bound (default: 5000)",false,5000,"T1 upper bound",cmd);
    TCLAP::ValueArg<double> upperBoundT2Arg("","upper-bound-t2","T2 value upper bound (default: 1000)",false,1000,"T2 upper bound",cmd);
    TCLAP::ValueArg<double> upperBoundM0Arg("","upper-bound-m0","M0 value upper bound (default: 5000)",false,5000,"M0 upper bound",cmd);
    TCLAP::ValueArg<double> lowerBoundB1Arg("","lower-bound-b1","B1 value lower bound (default: 0.2)",false,0.2,"B1 lower bound",cmd);
    TCLAP::ValueArg<double> upperBoundB1Arg("","upper-bound-b1","B1 value upper bound (default: 2)",false,2,"B1 upper bound",cmd);

    TCLAP::SwitchArg kMestimateSigmaArg("k","k-mest","M-estimation for K factor (experimental and apparently buggy)",cmd, false);
    TCLAP::ValueArg<double> b1SmoothingSigmaArg("s","b1-sigma","B1 smoothing sigma (in pixels, default: 3)",false,3,"B1 smoothing sigma",cmd);
    TCLAP::ValueArg<double> backgroundSignalThresholdArg("t","signal-thr","Background signal threshold (default: 100)",false,100,"Background signal threshold",cmd);

    TCLAP::ValueArg<unsigned int> numIterArg("n","num-iter","Number of optimization loops (default: 100)",false,100,"Number of iteration loops",cmd);
    TCLAP::ValueArg<unsigned int> numOptimizerIterArg("","opt-iter","Maximal number of optimizer iterations (default: 200)",false,200,"Maximal number of optimizer iterations",cmd);
    TCLAP::ValueArg<double> optimizerStopConditionArg("","opt-stop","Optimizer stopping threshold (default: 1.0e-4)",false,1.0e-4,"Optimizer stopping threshold",cmd);

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

    typedef itk::Image <double,3> InputImageType;
    typedef itk::Image <double,4> Image4DType;
    typedef InputImageType OutputImageType;
    typedef anima::CombinedRelaxometryEstimationImageFilter <InputImageType,OutputImageType> FilterType;

    FilterType::Pointer mainFilter = FilterType::New();

    // Read and set T1 relaxometry images
    bool readingFail = false;
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(t1Arg.getValue().c_str(),
                                                                           itk::IOFileModeEnum::ReadMode);

    if (!imageIO)
        readingFail = true;

    if (readingFail)
    {
        std::ifstream fileIn(t1Arg.getValue().c_str());
        if (!fileIn.is_open())
        {
            std::cerr << "Unable to read file: " << t1Arg.getValue() << std::endl;
            return EXIT_FAILURE;
        }

        while (!fileIn.eof())
        {
            char tmpStr[2048];
            fileIn.getline(tmpStr,2048);

            if (strcmp(tmpStr,"") == 0)
                continue;

            mainFilter->AddT1RelaxometryInput(anima::readImage <InputImageType> (tmpStr));
        }

        fileIn.close();
    }
    else
    {
        imageIO->SetFileName(t1Arg.getValue());
        imageIO->ReadImageInformation();

        // Do we have a (N+1)D image ? If not, exiting, otherwise considering the last dimension as each input
        unsigned int ndim = imageIO->GetNumberOfDimensions();

        if (ndim != (InputImageType::ImageDimension + 1))
        {
            std::cerr << "Unable to read file: " << t1Arg.getValue() << std::endl;
            return EXIT_FAILURE;
        }

        std::vector <InputImageType::Pointer> inputData;
        inputData = anima::getImagesFromHigherDimensionImage<Image4DType,InputImageType>(anima::readImage <Image4DType> (t1Arg.getValue()));

        for (unsigned int i = 0;i < inputData.size();++i)
            mainFilter->AddT1RelaxometryInput(inputData[i]);
    }

    std::ifstream inputT1Flip(t1FlipAnglesArg.getValue().c_str());
    if (!inputT1Flip.is_open())
    {
        std::cerr << "Text file for T1 flip angles not readable " << t1FlipAnglesArg.getValue() << std::endl;
        return EXIT_FAILURE;
    }

    std::vector <double> t1FlipAngles;
    while (!inputT1Flip.eof())
    {
        char tmpStr[8192];
        inputT1Flip.getline(tmpStr,8192);
        std::string workStr(tmpStr);
        if (workStr == "")
            continue;

        workStr.erase(workStr.find_last_not_of(" \n\r\t")+1);
        std::istringstream iss(workStr);
        std::string shortStr;
        iss >> shortStr;
        t1FlipAngles.push_back(std::stod(shortStr) * M_PI / 180.0);
    }

    mainFilter->SetT1FlipAngles(t1FlipAngles);
    mainFilter->SetTRT1Value(trT1Arg.getValue());

    // Read and set T2 relaxometry images
    unsigned int numT2Inputs = anima::setMultipleImageFilterInputsFromFileName<InputImageType,FilterType>(t2Arg.getValue(),mainFilter);

    mainFilter->SetT2FlipAngles(t2FlipAngleArg.getValue() * M_PI / 180.0,numT2Inputs);
    mainFilter->SetT2ExcitationFlipAngle(excitationT2FlipAngleArg.getValue() * M_PI / 180.0);
    mainFilter->SetT2EchoSpacing(t2EchoSpacingArg.getValue());
    mainFilter->SetTRT2Value(trT2Arg.getValue());
    mainFilter->SetKFactorMEstimation(kMestimateSigmaArg.isSet());

    mainFilter->SetM0UpperBound(upperBoundM0Arg.getValue());
    mainFilter->SetT1UpperBound(upperBoundT1Arg.getValue());
    mainFilter->SetT2UpperBound(upperBoundT2Arg.getValue());
    mainFilter->SetB1LowerBound(lowerBoundB1Arg.getValue());
    mainFilter->SetB1UpperBound(upperBoundB1Arg.getValue());

    if (maskArg.getValue() != "")
    {
        typedef itk::ImageFileReader < itk::Image <unsigned char, 3> > itkMaskReader;
        itkMaskReader::Pointer maskRead = itkMaskReader::New();
        maskRead->SetFileName(maskArg.getValue().c_str());
        maskRead->Update();

        mainFilter->SetComputationMask(maskRead->GetOutput());
    }

    mainFilter->SetAverageSignalThreshold(backgroundSignalThresholdArg.getValue());
    mainFilter->SetB1SmoothingSigma(b1SmoothingSigmaArg.getValue());

    mainFilter->SetNumberOfIterations(numIterArg.getValue());
    mainFilter->SetMaximumOptimizerIterations(numOptimizerIterArg.getValue());
    mainFilter->SetOptimizerStopCondition(optimizerStopConditionArg.getValue());

    mainFilter->SetNumberOfWorkUnits(nbpArg.getValue());

    itk::TimeProbe tmpTime;
    tmpTime.Start();

    mainFilter->Update();

    tmpTime.Stop();

    std::cout << "Total computation time: " << tmpTime.GetTotal() << std::endl;

    anima::writeImage <OutputImageType> (resT1Arg.getValue(),mainFilter->GetOutput(0));

    if (resT2Arg.getValue() != "")
        anima::writeImage <OutputImageType> (resT2Arg.getValue(),mainFilter->GetOutput(1));

    if (resT2M0Arg.getValue() != "")
        anima::writeImage <OutputImageType> (resT2M0Arg.getValue(),mainFilter->GetOutput(2));

    if (resT1M0Arg.getValue() != "")
    {
        double kFactor = mainFilter->GetKFactorM0();

        itk::ImageRegionIterator <OutputImageType> outputM0Iterator (mainFilter->GetOutput(2),mainFilter->GetOutput(2)->GetLargestPossibleRegion());

        while (!outputM0Iterator.IsAtEnd())
        {
            outputM0Iterator.Set(outputM0Iterator.Get() * kFactor);
            ++outputM0Iterator;
        }

        anima::writeImage <OutputImageType> (resT1M0Arg.getValue(),mainFilter->GetOutput(2));
    }

    if (resB1Arg.getValue() != "")
        anima::writeImage <OutputImageType> (resB1Arg.getValue(),mainFilter->GetOutput(3));

    if (resB1AddArg.getValue() != "")
        anima::writeImage <OutputImageType> (resB1AddArg.getValue(),mainFilter->GetOutput(4));

    std::cout << "Computed k factor between T2 and T1 M0 values: " << mainFilter->GetKFactorM0() << std::endl;

    return EXIT_SUCCESS;
}
