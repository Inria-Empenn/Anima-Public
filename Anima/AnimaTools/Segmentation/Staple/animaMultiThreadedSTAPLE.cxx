#include <iostream>
#include <tclap/CmdLine.h>

#include <animaMultiThreadedSTAPLEImageFilter.h>
#include <animaReadWriteFunctions.h>
#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkTimeProbe.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> lstArg("l","listname","File containing a list of segmentations",true,"","image list",cmd);
    TCLAP::ValueArg<std::string> resArg("o","outputname","Result image",true,"","result image",cmd);
    TCLAP::ValueArg<std::string> gtPriorArg("g","gtprior","Vector image of ground truth prior probabilities (for advanced use)",false,"","GT prior probabilities",cmd);
    TCLAP::ValueArg<std::string> presArg("s","outputparameters","Outputs expert parameters as an image",false,"","output parameters image",cmd);
    TCLAP::ValueArg<std::string> lmapArg("O","outputlabelmap","Outputs the classification as a label map",false,"","output label map",cmd);

    TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask",false,"","computation mask",cmd);
    TCLAP::ValueArg<unsigned int> maskDilRadiusArg("","mask-dil","Dilation radius to compute default work mask as dilated union of inputs, 0 means use the full image (default: 2)",false,2,"mask dilation radius",cmd);

    TCLAP::ValueArg<double> iniArg("d","initialparametervalue","Initial value of the diagonal parameters in STAPLE (default : 0.99)",false,0.99,"initial diagonal parameters",cmd);
    TCLAP::ValueArg<double> thrArg("t","relativethreshold","Relative stopping criterion for Staple (default : 1.0e-8)",false,1.0e-8,"stopping criterion",cmd);

    TCLAP::SwitchArg dmapArg("D","diagonalmap","Activate diagonal MAP estimator",cmd,false);
    TCLAP::SwitchArg fmapArg("F","fullmap","Activate full MAP estimator (overrides diagonal MAP)",cmd,false);

    TCLAP::SwitchArg mStrArg("M","missingstr","Account for missing structures. Used only if together with full MAP estimator",cmd,false);

    TCLAP::ValueArg<double> weightArg("w","mapweighting","Parameter to increase the weight of the prior in MAP estimation (default : 0.5)",false,0.5,"MAP weighting parameter",cmd);

    TCLAP::ValueArg<double> ampArg("a","alphadiagonalmap","Alpha parameter for the diagonal beta distributions (default : 5.0)",false,5.0,"diagonal alpha parameter",cmd);
    TCLAP::ValueArg<double> bmpArg("b","betadiagonalmap","Beta parameter for the diagonal beta distributions (default : 1.5)",false,1.5,"diagonal beta parameter",cmd);

    TCLAP::ValueArg<double> ampoArg("A","alphaoffdiagonalmap","Alpha parameter for the off-diagonal beta distributions (default : 1.5)",false,1.5,"off-diagonal alpha parameter",cmd);
    TCLAP::ValueArg<double> bmpoArg("B","betaoffdiagonalmap","Beta parameter for the off-diagonal beta distributions (default : 5.0)",false,5.0,"off-diagonal beta parameter",cmd);

    TCLAP::ValueArg<unsigned int> iterArg("i","iterations","Maximum number of iterations (default: 100)",false,100,"number of iterations",cmd);
    TCLAP::ValueArg<unsigned int> nbpArg("T","numberofthreads","Number of threads to run on (default : all available cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    int itMax = iterArg.getValue();
    double iniDiag = iniArg.getValue();
    double relThr = thrArg.getValue();

    unsigned int nbProcs = nbpArg.getValue();

    typedef itk::Image <unsigned short, 3> InputImageType;
    typedef itk::Image <double, 3> ParamsImageType;
    typedef anima::MultiThreadedSTAPLEImageFilter <InputImageType> StapleFilterType;
    typedef itk::VectorImage<double,3> OutputImageType;
    itk::TimeProbe tmpTime;
    tmpTime.Start();

    StapleFilterType::Pointer stapleFilter = StapleFilterType::New();
    stapleFilter->SetVerbose(true);
    stapleFilter->SetMAPWeighting(weightArg.getValue());
    stapleFilter->SetAccountForMissingStructures(mStrArg.getValue());
    stapleFilter->SetMaximumIterations(itMax);
    stapleFilter->SetRelativeConvergenceThreshold(relThr);

    if (maskArg.getValue() != "")
        stapleFilter->SetComputationMask(anima::readImage < itk::Image <unsigned char, 3> >(maskArg.getValue()));

    stapleFilter->SetMaskDilationRadius(maskDilRadiusArg.getValue());

    if (gtPriorArg.getValue() != "")
        stapleFilter->SetGTPriorImage(anima::readImage < itk::VectorImage <double,3> >(gtPriorArg.getValue()));

    std::ifstream fileIn(lstArg.getValue().c_str());
    if (!fileIn.is_open())
    {
        std::cerr << "Input file list " << lstArg.getValue() << " does not exist... Exiting..." << std::endl;
        return EXIT_FAILURE;
    }

    int nbPats = 0;

    while (!fileIn.eof())
    {
        char tmpStr[2048];
        fileIn.getline(tmpStr,2048);

        if (strcmp(tmpStr,"") == 0)
            continue;

        std::cout << "Loading image " << nbPats << " " << tmpStr << "..." << std::endl;
        stapleFilter->SetInput(nbPats,anima::readImage <InputImageType> (tmpStr));

        nbPats++;
    }

    fileIn.close();

    stapleFilter->InitializeNbClassesFromData();
    stapleFilter->InitializeExpertParameters(iniDiag);
    stapleFilter->InitializePriorFromData();

    if ((dmapArg.isSet()) || (fmapArg.isSet()))
    {
        double alphaMAP = ampArg.getValue();
        double betaMAP = bmpArg.getValue();
        double alphaMAPNonDiag = ampoArg.getValue();
        double betaMAPNonDiag = bmpoArg.getValue();

        if (dmapArg.isSet())
        {
            stapleFilter->SetMAPUpdate(StapleFilterType::DIAGONAL_MAP);
            stapleFilter->SetAlphaMAP(alphaMAP);
            stapleFilter->SetBetaMAP(betaMAP);
        }
        else
        {
            stapleFilter->SetMAPUpdate(StapleFilterType::FULL_MAP);
            stapleFilter->SetAlphaMAP(alphaMAP);
            stapleFilter->SetBetaMAP(betaMAP);
            stapleFilter->SetAlphaMAPNonDiag(alphaMAPNonDiag);
            stapleFilter->SetBetaMAPNonDiag(betaMAPNonDiag);
        }
    }
    else
        stapleFilter->SetMAPUpdate(StapleFilterType::STANDARD);

    stapleFilter->SetNumberOfWorkUnits(nbProcs);
    stapleFilter->Update();

    if (resArg.getValue() != "")
        anima::writeImage <OutputImageType> (resArg.getValue(),stapleFilter->GetOutput());

    if (lmapArg.getValue() != "")
        anima::writeImage <InputImageType> (lmapArg.getValue(),stapleFilter->GetClassificationAsLabelMap());

    if (presArg.getValue() == "")
        stapleFilter->PrintPerformanceParameters();
    else
        anima::writeImage <ParamsImageType> (presArg.getValue(),stapleFilter->GetExpertParametersAsImage());

    tmpTime.Stop();

    std::cout << "Total computation time: " << tmpTime.GetTotal() << std::endl;

    return EXIT_SUCCESS;
}
