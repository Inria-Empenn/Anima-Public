#include <cmath>

#include <animaMCMFileReader.h>
#include <animaReadWriteFunctions.h>

#include <itkTimeProbe.h>

#include <itkImage.h>
#include <animaMCMImage.h>

#include <itkImageRegionIterator.h>

#include <animaHyperbolicFunctions.h>

#include <tclap/CmdLine.h>

enum StringValue
{
    evNotDefined,
    evStringValue1,
    evStringValue2,
    evStringValue3,
    evStringValue4,
    evStringValue5,
    evStringValue6,
    evStringValue7,
    evStringValue8,
    evStringValue9,
    evStringValue10,
    evStringValue11
};

void InitializeStringArguments(std::map<std::string, StringValue> &mapStringValues)
{
    mapStringValues["fa"] = evStringValue1;
    mapStringValues["md"] = evStringValue2;
    mapStringValues["dpara"] = evStringValue3;
    mapStringValues["dperp"] = evStringValue4;
    mapStringValues["vr"] = evStringValue5;
    mapStringValues["cs"] = evStringValue6;
    mapStringValues["kpara"] = evStringValue7;
    mapStringValues["kperp"] = evStringValue8;
    mapStringValues["od"] = evStringValue9;
    mapStringValues["eaf"] = evStringValue10;
    mapStringValues["wiso"] = evStringValue11;

    // cl is equivalent to fa under cylindrical symmetry
    // cp is zero under cylindrical symmetry
}

int main(int argc,  char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","ddi","DDI volume",true,"","DDI volume",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","Result scalar image",true,"","result scalar image",cmd);

    TCLAP::SwitchArg averageArg("A","average-val","Compute map at voxel level",cmd,false);

    TCLAP::ValueArg<std::string> typeArg("t","scalartype","Please choose between\nFractional Anisotropy (FA)\nMean Diffusivity (MD)\nAxial diffusivity (DPARA)\nRadial Diffusivity (DPERP)\nVolume Ratio (VR)\nSpherical Coefficient (CS) \nExtraAxionalFraction (EAF) \n Free water (wiso)" ,true,"","scalar type",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    std::map<std::string, StringValue> mapStringValues;
    InitializeStringArguments(mapStringValues);

    typedef double ScalarType;
    typedef anima::MCMImage <ScalarType, 3> TInputImageType;
    typedef itk::Image <ScalarType,3> TOutputImageType;

    /** Image typedef support */
    typedef TOutputImageType::Pointer OutputImagePointer;

    anima::MCMFileReader <double,3> mcmReader;
    mcmReader.SetFileName(inArg.getValue());
    mcmReader.Update();

    typedef itk::ImageRegionIterator <TInputImageType> InputIteratorType;
    typedef itk::ImageRegionIterator <TOutputImageType> OutputIteratorType;

    TInputImageType::Pointer inputImage = mcmReader.GetModelVectorImage();
    anima::MultiCompartmentModel::Pointer mcm = inputImage->GetDescriptionModel();

    OutputImagePointer outImage = TOutputImageType::New();
    outImage->Initialize();
    outImage->SetRegions(inputImage->GetLargestPossibleRegion());
    outImage->SetOrigin(inputImage->GetOrigin());
    outImage->SetDirection(inputImage->GetDirection());
    outImage->SetSpacing(inputImage->GetSpacing());
    outImage->Allocate();

    bool hasFreeWater = (mcm->GetCompartment(0)->GetCompartmentType() == anima::FreeWater);
    if (!hasFreeWater)
    {
        std::cerr << "We handle only models with a free water component for now" << std::endl;
        return -1;
    }

    unsigned int numberOfCompartments = mcm->GetNumberOfCompartments() - 1;

    InputIteratorType ddiItr (inputImage,inputImage->GetLargestPossibleRegion());
    OutputIteratorType outItr (outImage,outImage->GetLargestPossibleRegion());

    std::transform(typeArg.getValue().begin(),typeArg.getValue().end(),typeArg.getValue().begin(), ::tolower);

    while (!outItr.IsAtEnd())
    {
        mcm->SetModelVector(ddiItr.Get());
        double wIso = mcm->GetCompartmentWeight(0);
        double dIso = mcm->GetCompartment(0)->GetAxialDiffusivity();

        double voxelAxialDiffusivity = 0;
        double voxelRadialDiffusivity = 0;
        double voxelAxialKurtosis = 0;
        double voxelRadialKurtosis = 0;
        double meanKappa = 0;
        double extraAxonalFraction = 0;

        for (unsigned int m = 0;m < numberOfCompartments;++m)
        {
            double kappaVal = mcm->GetCompartment(m+1)->GetOrientationConcentration();
            double xiVal = anima::xi(kappaVal);

            double dVal = mcm->GetCompartment(m+1)->GetAxialDiffusivity();
            double aVal = mcm->GetCompartment(m+1)->GetExtraAxonalFraction();
            double wVal = mcm->GetCompartmentWeight(m+1);

            voxelAxialDiffusivity += wVal * dVal * (1.0 - 2.0 * aVal * xiVal);
            voxelRadialDiffusivity += wVal * dVal * ((1.0 - aVal) / (kappaVal + 1.0) + aVal * xiVal);

            voxelAxialKurtosis += wVal * dVal * dVal * (3.0 * (1.0 - aVal) * (1.0 - aVal)
                                                        + 6.0 * aVal * (1.0 - aVal) * anima::jtwo(kappaVal)
                                                        + aVal * aVal * anima::jfour(kappaVal));
            voxelRadialKurtosis += wVal * dVal * dVal * (aVal * aVal + 8.0 * aVal * (1.0 - aVal) / (kappaVal + 1.0) + 8.0 * (1.0 - aVal) * (1.0 - aVal) / ((kappaVal + 1.0) * (kappaVal + 1.0))
                                                         - 2.0 * (aVal * aVal + 4.0 * aVal * (1.0 - aVal) / (kappaVal + 1.0)) * anima::jtwo(kappaVal)
                                                         + aVal * aVal * anima::jfour(kappaVal));

            double thrKappa = (wIso > 0.99) ? 0 : kappaVal;
            meanKappa += wVal * thrKappa;

            extraAxonalFraction += wVal * aVal;
        }
        
        if (extraAxonalFraction < 1.0e-6)
            extraAxonalFraction = 0;
        else
            extraAxonalFraction /= (1.0 - wIso);

        if (averageArg.isSet())
        {
            voxelAxialDiffusivity += wIso * dIso;
            voxelRadialDiffusivity += wIso * dIso;

            voxelAxialKurtosis += 3.0 * wIso * dIso * dIso;
            voxelRadialKurtosis += 3.0 * wIso * dIso * dIso;
        }
        else
        {
            voxelAxialDiffusivity /= (1.0 - wIso);
            voxelRadialDiffusivity /= (1.0 - wIso);

            voxelAxialKurtosis /= (1.0 - wIso);
            voxelRadialKurtosis /= (1.0 - wIso);

            meanKappa /= (1.0 - wIso);
        }

        double scalarVal = 0;
        double diffVal = voxelAxialDiffusivity - voxelRadialDiffusivity;
        double normVal = std::sqrt(voxelAxialDiffusivity * voxelAxialDiffusivity + 2.0 * voxelRadialDiffusivity * voxelRadialDiffusivity);

        switch (mapStringValues[typeArg.getValue()])
        {
            case evStringValue1:
                scalarVal = diffVal / normVal;
                break;

            case evStringValue2:
                scalarVal = (voxelAxialDiffusivity + 2.0 * voxelRadialDiffusivity) / 3.0;
                break;

            case evStringValue3:
                scalarVal = voxelAxialDiffusivity;
                break;

            case evStringValue4:
                scalarVal = voxelRadialDiffusivity;
                break;

            case evStringValue5:
                scalarVal = std::sqrt(voxelAxialDiffusivity) * voxelRadialDiffusivity / std::pow(dIso, 1.5);
                break;

            case evStringValue6:
                scalarVal = 1.5 * voxelRadialDiffusivity / normVal;
                break;

            case evStringValue7:
                scalarVal = voxelAxialKurtosis / (voxelAxialDiffusivity * voxelAxialDiffusivity) - 3.0;
                break;

            case evStringValue8:
                scalarVal = voxelRadialKurtosis / (voxelRadialDiffusivity * voxelRadialDiffusivity) - 3.0;
                break;

            case evStringValue9:
                scalarVal = 2.0 / M_PI * std::atan(1.0 / meanKappa);
                break;

            case evStringValue10:
                scalarVal = extraAxonalFraction;
                break;

            case evStringValue11:
                scalarVal = wIso;
            break;

            default:
                std::cerr << "Unrecognized scalar type map. Please choose between\nFractional Anisotropy (FA)\nMean Diffusivity (MD)\nAxial diffusivity (DPARA)\nRadial Diffusivity (DPERP)\nVolume Ratio (VR)\nSpherical Coefficient (CS)" << std::endl;
                exit(-1);
        }

        outItr.Set(scalarVal);

        ++ddiItr;
        ++outItr;
    }

    anima::writeImage <TOutputImageType> (outArg.getValue(),outImage);

    return 0;
}
