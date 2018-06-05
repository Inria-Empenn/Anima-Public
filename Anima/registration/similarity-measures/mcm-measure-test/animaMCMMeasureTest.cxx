#include <animaMultiCompartmentModel.h>
#include <animaMCMPairingMeanSquaresImageToImageMetric.h>
#include <animaMCMMeanSquaresImageToImageMetric.h>
#include <animaMCMCorrelationImageToImageMetric.h>
#include <animaMTPairingCorrelationImageToImageMetric.h>
#include <itkImageRegionIterator.h>
#include <animaLogRigid3DTransform.h>
#include <itkTimeProbe.h>

#include <fstream>
#include <tclap/CmdLine.h>
#include <animaGradientFileReader.h>
#include <animaMCMLinearInterpolateImageFunction.h>
#include <animaMCMFileReader.h>

int main(int ac, const char** av)
{
    // Parsing arguments
    TCLAP::CmdLine  cmd("INRIA / IRISA - VisAGeS Team", ' ', ANIMA_VERSION);

    // Setting up parameters
    TCLAP::ValueArg<std::string> inArg("i","im","input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> bvalArg("b","bvals","B-values",true,"","b-values",cmd);
    TCLAP::ValueArg<std::string> bvecArg("v","bvec","Gradient direction",true,"","gradient directions",cmd);
    TCLAP::ValueArg<std::string> outArg("o","out","Output file",true,"","output file",cmd);

    try
    {
        cmd.parse(ac,av);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef anima::MCMFileReader <double,3> MCMReaderType;
    MCMReaderType refReader;
    refReader.SetFileName(inArg.getValue());
    refReader.Update();

    MCMReaderType movingReader;
    movingReader.SetFileName(inArg.getValue());
    movingReader.Update();

    typedef anima::GradientFileReader <vnl_vector_fixed <double,3>, double> GradientReaderType;
    GradientReaderType gradientReader;
    gradientReader.SetGradientFileName(bvecArg.getValue());
    gradientReader.SetBValueBaseString(bvalArg.getValue());
    gradientReader.SetGradientIndependentNormalization(true);

    gradientReader.Update();

    typedef anima::MCMPairingMeanSquaresImageToImageMetric <double, double, 3> BasicMetricType;
    typedef BasicMetricType::Pointer BasicMetricPointer;

    BasicMetricPointer testBasicMetric = BasicMetricType::New();
    testBasicMetric->SetFixedImage(refReader.GetModelVectorImage());
    testBasicMetric->SetMovingImage(movingReader.GetModelVectorImage());

    typedef MCMReaderType::OutputImageType ImageType;
    ImageType::RegionType region;
    region.SetIndex(0,68);
    region.SetIndex(1,56);
    region.SetIndex(2,62);
    region.SetSize(0,5);
    region.SetSize(1,5);
    region.SetSize(2,5);
    testBasicMetric->SetFixedImageRegion(region);

    typedef anima::MCMLinearInterpolateImageFunction <ImageType> InterpolateFunctionType;
    InterpolateFunctionType::Pointer interpolator = InterpolateFunctionType::New();
    interpolator->SetInputImage(movingReader.GetModelVectorImage());
    interpolator->SetReferenceOutputModel(refReader.GetModelVectorImage()->GetDescriptionModel());
    testBasicMetric->SetInterpolator(interpolator);

    typedef anima::MCMMeanSquaresImageToImageMetric <double, double, 3> MCMMetricType;
    typedef MCMMetricType::Pointer MCMMetricPointer;
    MCMMetricPointer testMetric = MCMMetricType::New();

    testMetric->SetFixedImage(refReader.GetModelVectorImage());
    testMetric->SetMovingImage(movingReader.GetModelVectorImage());
    testMetric->SetFixedImageRegion(region);
    // Use default small and large delta values
    testMetric->SetGradientStrengths(gradientReader.GetGradientStrengths());
    testMetric->SetGradientDirections(gradientReader.GetGradients());
    testMetric->SetInterpolator(interpolator);

    typedef anima::MCMCorrelationImageToImageMetric <double, double, 3> MCMCorrelationMetricType;
    typedef MCMCorrelationMetricType::Pointer MCMCorrelationMetricPointer;
    MCMCorrelationMetricPointer testCorrelationMetric = MCMCorrelationMetricType::New();

    testCorrelationMetric->SetFixedImage(refReader.GetModelVectorImage());
    testCorrelationMetric->SetMovingImage(movingReader.GetModelVectorImage());
    testCorrelationMetric->SetFixedImageRegion(region);
    // Use default small and large delta values
    testCorrelationMetric->SetGradientStrengths(gradientReader.GetGradientStrengths());
    testCorrelationMetric->SetGradientDirections(gradientReader.GetGradients());
    testCorrelationMetric->SetInterpolator(interpolator);

    typedef anima::LogRigid3DTransform <double> TransformType;
    TransformType::Pointer trsf = TransformType::New();
    trsf->SetIdentity();
    ImageType::IndexType centerIndex;
    centerIndex[0] = 70;
    centerIndex[1] = 58;
    centerIndex[2] = 64;
    TransformType::InputPointType center;
    movingReader.GetModelVectorImage()->TransformIndexToPhysicalPoint(centerIndex,center);
    trsf->SetCenter(center);
    testBasicMetric->SetTransform(trsf);
    testBasicMetric->SetModelRotation(BasicMetricType::Superclass::PPD);
    testBasicMetric->PreComputeFixedValues();
    testMetric->SetTransform(trsf);
    testMetric->SetModelRotation(MCMMetricType::Superclass::PPD);
    testMetric->PreComputeFixedValues();
    testCorrelationMetric->SetTransform(trsf);
    testCorrelationMetric->SetModelRotation(MCMMetricType::Superclass::PPD);
    testCorrelationMetric->PreComputeFixedValues();

    TransformType::ParametersType parameters = trsf->GetParameters();
    std::ofstream outFile(outArg.getValue().c_str());

    itk::TimeProbe timerL2Approx, timerL2, timerPairingOneToOne, timerPairing, timerCorrelation, timerCorrelationApprox;
    outFile.precision(10);
    //for (int i = -90;i <= 90;++i)
    for (int i = -100;i <= 100;++i)
    {
        //double angle = M_PI * i / 180.0;
        //parameters[2] = angle;
        double translation = i / 10.0;
        parameters[4] = translation;

        timerL2.Start();
        testMetric->SetForceApproximation(false);
        double metricValue = testMetric->GetValue(parameters);
        timerL2.Stop();

        timerL2Approx.Start();
        testMetric->SetForceApproximation(true);
        double metricValueApprox = testMetric->GetValue(parameters);
        timerL2Approx.Stop();

        timerCorrelation.Start();
        testCorrelationMetric->SetForceApproximation(false);
        double metricCorrelationValue = testCorrelationMetric->GetValue(parameters);
        timerCorrelation.Stop();

        timerCorrelationApprox.Start();
        testCorrelationMetric->SetForceApproximation(true);
        double metricApproxCorrelationValue = testCorrelationMetric->GetValue(parameters);
        timerCorrelationApprox.Stop();

        timerPairingOneToOne.Start();
        testBasicMetric->SetOneToOneMapping(true);
        double basicMetricValueOneToOne = testBasicMetric->GetValue(parameters);
        timerPairingOneToOne.Stop();

        timerPairing.Start();
        testBasicMetric->SetOneToOneMapping(false);
        double basicMetricValue = testBasicMetric->GetValue(parameters);
        timerPairing.Stop();

        outFile << translation << " " << metricValue << " " << metricValueApprox << " "
                << basicMetricValue << " " << basicMetricValueOneToOne << " "
                << metricCorrelationValue << " " << metricApproxCorrelationValue << std::endl;
    }

    std::cout << "Time L2: " << timerL2.GetTotal() << std::endl;
    std::cout << "Time L2 approx: " << timerL2Approx.GetTotal() << std::endl;
    std::cout << "Time pairing: " << timerPairing.GetTotal() << std::endl;
    std::cout << "Time pairing one to one: " << timerPairingOneToOne.GetTotal() << std::endl;
    std::cout << "Time correlation: " << timerCorrelation.GetTotal() << std::endl;
    std::cout << "Time correlation approx: " << timerCorrelationApprox.GetTotal() << std::endl;

    outFile.close();

    return EXIT_SUCCESS;
}
