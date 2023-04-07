#include <tclap/CmdLine.h>

#include <animaReadWriteFunctions.h>

#include <animaVectorOperations.h>
#include <animaPrivateMultiCompartmentModelCreator.h>
#include <animaMCMFileWriter.h>

using namespace anima;

int main(int argc, char *argv[])
{
    // Parsing arguments
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    // Setting up parameters
    TCLAP::ValueArg<std::string> inArg("i","input","list of text file",true,"","text file with value of a fascicule",cmd);
    TCLAP::ValueArg<std::string> resArg("o","output","grid of these different value",true,"","grid of these different value",cmd);
    TCLAP::ValueArg<unsigned int> hArg("", "scale-0", "number of horizontal value", false, 11, "scalar type", cmd);
    TCLAP::ValueArg<unsigned int> vArg("", "scale-1", "number of vertical value", false, 11, "scalar", cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        exit(-1);
    }

    std::ifstream inputFile(inArg.getValue().c_str());

    if (!inputFile.is_open())
    {
        std::cerr << "Please provide usable file with input DDIs" << std::endl;
        abort();
    }

    unsigned int numInput = 0;
    unsigned int nbPoints0 = hArg.getValue();
    unsigned int nbPoints1 = vArg.getValue();

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

    anima::PrivateMultiCompartmentModelCreator mcmCreator;
    mcmCreator.SetCompartmentType(anima::DDI);
    mcmCreator.SetNumberOfCompartments(numInput);
    mcmCreator.SetModelWithFreeWaterComponent(false);

    anima::MultiCompartmentModel::Pointer mcm = mcmCreator.GetNewMultiCompartmentModel();

    ImageType::RegionType region(start, size);
    ImageType::Pointer image = ImageType::New();
    image->SetRegions(region);
    image->SetVectorLength(mcm->GetSize());
    image->Allocate();

    ImageType::IndexType pixelIndex;
    pixelIndex[2] = 0;

    std::vector<ScalarType> weight(4);
    ScalarType scale0 = 1.0/(nbPoints0 - 1.0);
    ScalarType scale1 = 1.0/(nbPoints1 - 1.0);

    ImageType::PixelType  voxelOutputValue(mcm->GetSize());
    voxelOutputValue.Fill(0);
    mcm->SetModelVector(voxelOutputValue);

    for (unsigned int k = 0; k < numInput; k++)
    {
        anima::BaseCompartment *workCompartment = mcm->GetCompartment(k);

        workCompartment->SetAxialDiffusivity(d_fascicule[k]);
        workCompartment->SetOrientationConcentration(kappa[k]);
        workCompartment->SetExtraAxonalFraction(v[k]);

        anima::TransformCartesianToSphericalCoordinates(directions[k],directions[k]);
        workCompartment->SetOrientationTheta(directions[k][0]);
        workCompartment->SetOrientationPhi(directions[k][1]);
    }

    for (unsigned int i = 0; i < nbPoints0; i++)
    {
        for (unsigned int j = 0; j < nbPoints1; j++)
        {       
            pixelIndex[0] = i;
            pixelIndex[1] = j;

            weight[0] = (1 - scale0*i) * (1 - scale1*j);
            weight[1] = (1 - scale0*i) * scale1*j;
            weight[2] = scale0*i * (1 - scale1*j);
            weight[3] = scale0 * i * scale1 * j;

            mcm->SetCompartmentWeights(weight);
            voxelOutputValue = mcm->GetModelVector();

            //if (((i == 0)&&((j == 0)||(j == nbPoints1 - 1))) || ((i == nbPoints0 - 1)&&((j == 0)||(j == nbPoints1 - 1))))
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
