#include <tclap/CmdLine.h>
#include <cmath>

#include <animaDDIAveragingTools.h>
#include <animaMCMFileWriter.h>
#include <animaMCMImage.h>

#include <animaMultiCompartmentModelCreator.h>

int main(int argc, char *argv[])
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","inputddi","list of text file",true,"","text file with value of a fascicule",cmd);
    TCLAP::ValueArg<std::string> resArg("o","output","extrapolate grid from input pixels",true,"","extrapolate grid from input pixels",cmd);
    TCLAP::ValueArg<double> hArg("", "scale-0", "number of horizontal value", false, 11, "scalar type", cmd);
    TCLAP::ValueArg<double> vArg("", "scale-1", "number of vertical value", false, 11, "scalar", cmd);
    TCLAP::ValueArg<int> methodArg("m", "method", "method use to average, Classic : 0, Tensor : 1, logVMF : 2, covarianceAnalytic : 3", true , 3, "method use to average", cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    std::ifstream inputFile(inArg.getValue().c_str());

    if (!inputFile.is_open())
    {
        std::cerr << "Please provide usable file with input DDIs" << std::endl;
        abort();
    }

    int method = methodArg.getValue();
    double  nbPoints0 = hArg.getValue();
    double  nbPoints1 = vArg.getValue();
    unsigned int numInput = 0;

    typedef double ScalarType;
    std::vector<vnl_vector<ScalarType> > directions;
    vnl_vector<ScalarType> onedirection(3);
    std::vector<ScalarType> kappa;
    std::vector<ScalarType> d_fascicule;
    std::vector<ScalarType> v;

    ScalarType temp0, temp1, temp2, temp3, temp4, temp5;
    while (!inputFile.eof())
    {
        char tmpStr[2048];
        inputFile.getline(tmpStr,2048);

        if (strcmp(tmpStr,"") == 0)
            continue;

        std::cout << "Loading pixel " << numInput << ": " << tmpStr << std::endl;

        std::ifstream pixelFile(tmpStr);

        if (!pixelFile.is_open())
        {
            std::cerr << "Please provide usable pixel file" << std::endl;
            abort();
        }
        pixelFile >> temp0 >> temp1 >> temp2 >> temp3 >> temp4 >> temp5;

        onedirection(0) = temp0;
        onedirection(1) = temp1;
        onedirection(2) = temp2;
        directions.push_back(onedirection);
        kappa.push_back(temp3);
        d_fascicule.push_back(temp4);
        v.push_back(temp5);

        numInput++;
    }

    typedef anima::MCMImage <ScalarType, 3> ImageType;
    ImageType::IndexType start;
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;

    ImageType::SizeType size;
    size[0] = nbPoints0;
    size[1] = nbPoints1;
    size[2] = 1;

    anima::MultiCompartmentModelCreator mcmCreator;
    mcmCreator.SetCompartmentType(anima::DDI);
    mcmCreator.SetNumberOfCompartments(1);
    mcmCreator.SetModelWithFreeWaterComponent(false);

    anima::MultiCompartmentModel::Pointer mcm = mcmCreator.GetNewMultiCompartmentModel();

    ImageType::RegionType region(start, size);
    ImageType::Pointer image = ImageType::New();
    image->SetRegions(region);
    image->SetVectorLength(mcm->GetSize());
    image->Allocate();

    ImageType::IndexType pixelIndex;
    pixelIndex[2] = 0;

    ScalarType scale0 = 1.0 / (nbPoints0 - 1);
    ScalarType scale1 = 1.0 / (nbPoints1 - 1);

    ImageType::PixelType voxelOutputValue(mcm->GetSize());
    std::vector <double> weights(4,1);
    std::vector <double> outputWeights(1,1);
    mcm->SetCompartmentWeights(outputWeights);

    for (int i = 0; i < nbPoints0; i++)
    {
        for (int j = 0; j < nbPoints1; j++)
        {
            pixelIndex[0] = i;
            pixelIndex[1] = j;
            weights[0] = (1 - scale0*i) * (1 - scale1*j);
            weights[1] = (1 - scale0*i) * scale1*j;
            weights[2] = scale0*i * (1 - scale1*j);
            weights[3] = scale0 * i * scale1 * j;

            double averageNu = 0, averageDiffusivity = 0, averageKappa = 0, averageWeight = 1;

            anima::DDIAveraging(v,d_fascicule,kappa,directions,weights,method,averageNu,averageDiffusivity,averageKappa,onedirection,averageWeight);

            anima::BaseCompartment *workCompartment = mcm->GetCompartment(0);
            workCompartment->SetAxialDiffusivity(averageDiffusivity);
            workCompartment->SetOrientationConcentration(averageKappa);
            workCompartment->SetExtraAxonalFraction(averageNu);
            anima::TransformCartesianToSphericalCoordinates(onedirection,onedirection);
            workCompartment->SetOrientationTheta(onedirection[0]);
            workCompartment->SetOrientationPhi(onedirection[1]);

            voxelOutputValue = mcm->GetModelVector();
            image->SetPixel(pixelIndex, voxelOutputValue);
        }
    }

    image->SetDescriptionModel(mcm);
    anima::MCMFileWriter <ScalarType,3> mcmWriter;
    mcmWriter.SetFileName(resArg.getValue());
    mcmWriter.SetInputImage(image);

    mcmWriter.Update();

    return 0;
}
