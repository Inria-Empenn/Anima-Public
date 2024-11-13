#include <animaGradientFileReader.h>
#include <animaODFEstimatorImageFilter.h>
#include <animaReadWriteFunctions.h>

#include <itkTimeProbe.h>

#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ', ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg(
        "i", "input",
        "List of diffusion weighted images or 4D volume",
        true, "", "input diffusion images", cmd);
    TCLAP::ValueArg<std::string> resArg(
        "o", "outputfile",
        "Result ODF image",
        true, "", "result ODF image", cmd);
    TCLAP::ValueArg<std::string> b0OutArg(
        "O", "output-b0",
        "output_b0",
        false, "", "result B0 image", cmd);
    TCLAP::ValueArg<std::string> varOutArg(
        "V", "output-variance",
        "output_variance",
        false, "", "result noise variance image", cmd);

    TCLAP::ValueArg<std::string> gradArg(
        "g", "gradientlist",
        "List of gradients (text file)",
        true, "", "list of gradients", cmd);
    TCLAP::ValueArg<std::string> bvalArg(
        "b", "bval",
        "input b-values (for checking single shell)",
        true, "", "Input b-values", cmd);
    TCLAP::ValueArg<std::string> refB0Arg(
        "r", "ref-b0",
        "Externally estimated B0 (otherwise mean of B0s is used)",
        false, "", "External B0 image", cmd);
    TCLAP::SwitchArg bvalueScaleArg(
        "B", "b-no-scale",
        "Do not scale b-values according to gradient norm",
        cmd);
    TCLAP::ValueArg<int> selectedBvalArg("v", "select-bval", "B-value shell used to estimate ODFs (default: first one in data volume above 10)", false, -1, "b-value shell selection", cmd);

    TCLAP::ValueArg<double> lambdaArg(
        "l", "lambda",
        "Lambda regularization parameter (see Descoteaux MRM 2007)",
        false, 0.006, "lambda for regularization", cmd);
    TCLAP::ValueArg<unsigned int> orderArg(
        "k", "order",
        "Order of spherical harmonics basis",
        false, 4, "Order of SH basis", cmd);

    TCLAP::ValueArg<double> sharpFactorArg(
        "s", "sharpenratio",
        "Ratio for sharpening ODFs (see Descoteaux TMI 2009, default : 0.255)",
        false, 0.255, "sharpening ratio", cmd);
    TCLAP::SwitchArg sharpenArg(
        "S", "sharpenodf",
        "Sharpen ODF ? (default: no)", cmd, false);

    TCLAP::ValueArg<std::string> normSphereArg(
        "n", "normalizefile",
        "Sphere tesselation for normalization",
        false, "", "Normalization sphere file", cmd);
    TCLAP::SwitchArg normalizeArg(
        "N", "normalize",
        "Normalize ODF ? (default: no)",
        cmd, false);

    TCLAP::SwitchArg radialArg(
        "R", "radialestimation",
        "Use radial estimation (see Aganj et al) ? (default: no)",
        cmd, false);
    TCLAP::ValueArg<double> aganjRegFactorArg(
        "d", "adr",
        "Delta threshold for signal regularization, only use if R option activated (see Aganj et al, default : 0.001)",
        false, 0.001, "delta signal regularization", cmd);

    TCLAP::ValueArg<unsigned int> nbpArg(
        "T", "nb-threads",
        "An integer value specifying the number of threads to run on (default: all cores).",
        false, itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(), "number of threads", cmd);

    try
    {
        cmd.parse(argc, argv);
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    using MainFilterType = anima::ODFEstimatorImageFilter<double, double>;
    using InputImageType = MainFilterType::TInputImage;

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

    anima::setMultipleImageFilterInputsFromFileName<InputImageType, MainFilterType>(inArg.getValue(), mainFilter);

    if (refB0Arg.getValue() != "")
        mainFilter->SetReferenceB0Image(anima::readImage<InputImageType>(refB0Arg.getValue()));

    using GFReaderType = anima::GradientFileReader<std::vector<double>, double>;
    GFReaderType gfReader;
    gfReader.SetGradientFileName(gradArg.getValue());
    gfReader.SetBValueBaseString(bvalArg.getValue());
    gfReader.SetGradientIndependentNormalization(bvalueScaleArg.isSet());
    gfReader.SetB0ValueThreshold(10);

    gfReader.Update();

    GFReaderType::GradientVectorType directions = gfReader.GetGradients();
    GFReaderType::BValueVectorType mb = gfReader.GetBValues();

    for (unsigned int i = 0; i < directions.size(); ++i)
        mainFilter->AddGradientDirection(i, directions[i]);

    mainFilter->SetBValuesList(mb);
    mainFilter->SetBValueShellSelected(selectedBvalArg.getValue());
    mainFilter->SetNumberOfWorkUnits(nbpArg.getValue());

    itk::TimeProbe tmpTime;
    tmpTime.Start();
    mainFilter->Update();
    tmpTime.Stop();

    std::cout << "\nExecution Time: " << tmpTime.GetTotal() << "s" << std::endl;

    anima::writeImage<MainFilterType::TOutputImage>(resArg.getValue(), mainFilter->GetOutput());

    if (b0OutArg.getValue() != "")
        anima::writeImage<InputImageType>(b0OutArg.getValue(), mainFilter->GetEstimatedB0Image());

    if (varOutArg.getValue() != "")
        anima::writeImage<InputImageType>(varOutArg.getValue(), mainFilter->GetEstimatedVarianceImage());

    return EXIT_SUCCESS;
}
