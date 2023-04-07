#include <tclap/CmdLine.h>
#include <animaMCMImage.h>
#include <animaMCML2DistanceComputer.h>
#include <animaBaseTensorTools.h>

#include <animaMCMPrivateFileReader.h>
#include <animaReadWriteFunctions.h>
#include <itkImageRegionIterator.h>

int main(int argc, char **argv)
{
    typedef anima::MCMImage <double,3> ImageType;
    typedef ImageType::MCMPointer MCMPointer;
    typedef anima::MCMPrivateFileReader <double,3> ImageReaderType;

    TCLAP::CmdLine cmd("Computes various squared distances between estimated and reference multi-tensors models.\n\
                       Assumes same number of compartments at each voxel.\n INRIA / IRISA - VisAGeS/Empenn Team",
                       ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> testArg("t","test","Test MCM data",true,"","test MCM",cmd);
    TCLAP::ValueArg<std::string> refArg("r","ref","Reference MCM data",true,"","reference MCM",cmd);

    TCLAP::ValueArg<std::string> outArg("o","out-dist","output squared tensor distance image",true,"","output squared tensor distance image",cmd);
    TCLAP::ValueArg<std::string> outWSEArg("O","out-wse","output squaredWSE image",true,"","output squared WSE image",cmd);
    TCLAP::ValueArg<std::string> outWSEIndArg("W","out-wse-ind","output squared WSE image (independent of tensor pairing)",true,"","output squared WSE independent image",cmd);
    TCLAP::ValueArg<std::string> outL2Arg("l","out-l2","output squared L2 norm image",true,"","output squared L2 norm image",cmd);

    TCLAP::ValueArg<std::string> maskArg("m","mask","mask image",true,"","mask image",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    typedef itk::Image <unsigned char, 3> MaskImageType;
    typedef itk::Image <double, 3> OutputImageType;
    typedef OutputImageType::Pointer OutputImagePointer;

    MaskImageType::Pointer maskImage = anima::readImage <MaskImageType> (maskArg.getValue());

    ImageReaderType refRead;
    refRead.SetFileName(refArg.getValue());
    refRead.Update();

    ImageType::Pointer refImage = refRead.GetModelVectorImage();
    MCMPointer refMCM = refImage->GetDescriptionModel()->Clone();

    ImageReaderType testRead;
    testRead.SetFileName(testArg.getValue());
    testRead.Update();

    ImageType::Pointer testImage = testRead.GetModelVectorImage();
    MCMPointer testMCM = testImage->GetDescriptionModel()->Clone();

    OutputImagePointer wseImage = OutputImageType::New();
    wseImage->Initialize();
    wseImage->SetOrigin(maskImage->GetOrigin());
    wseImage->SetSpacing(maskImage->GetSpacing());
    wseImage->SetDirection(maskImage->GetDirection());
    wseImage->SetRegions(maskImage->GetLargestPossibleRegion());
    wseImage->Allocate();
    wseImage->FillBuffer(0);

    OutputImagePointer wseIndImage = OutputImageType::New();
    wseIndImage->Initialize();
    wseIndImage->SetOrigin(maskImage->GetOrigin());
    wseIndImage->SetSpacing(maskImage->GetSpacing());
    wseIndImage->SetDirection(maskImage->GetDirection());
    wseIndImage->SetRegions(maskImage->GetLargestPossibleRegion());
    wseIndImage->Allocate();
    wseIndImage->FillBuffer(0);

    OutputImagePointer distImage = OutputImageType::New();
    distImage->Initialize();
    distImage->SetOrigin(maskImage->GetOrigin());
    distImage->SetSpacing(maskImage->GetSpacing());
    distImage->SetDirection(maskImage->GetDirection());
    distImage->SetRegions(maskImage->GetLargestPossibleRegion());
    distImage->Allocate();
    distImage->FillBuffer(0);

    OutputImagePointer l2Image = OutputImageType::New();
    l2Image->Initialize();
    l2Image->SetOrigin(maskImage->GetOrigin());
    l2Image->SetSpacing(maskImage->GetSpacing());
    l2Image->SetDirection(maskImage->GetDirection());
    l2Image->SetRegions(maskImage->GetLargestPossibleRegion());
    l2Image->Allocate();
    l2Image->FillBuffer(0);

    itk::ImageRegionIterator <MaskImageType> maskItr(maskImage, maskImage->GetLargestPossibleRegion());
    itk::ImageRegionIterator <OutputImageType> wseItr(wseImage, maskImage->GetLargestPossibleRegion());
    itk::ImageRegionIterator <OutputImageType> wseIndItr(wseIndImage, maskImage->GetLargestPossibleRegion());
    itk::ImageRegionIterator <OutputImageType> distItr(distImage, maskImage->GetLargestPossibleRegion());
    itk::ImageRegionIterator <OutputImageType> l2Itr(l2Image, maskImage->GetLargestPossibleRegion());

    itk::ImageRegionIterator <ImageType> testItr(testImage, maskImage->GetLargestPossibleRegion());
    itk::ImageRegionIterator <ImageType> refItr(refImage, maskImage->GetLargestPossibleRegion());

    unsigned int numCompartments = testMCM->GetNumberOfCompartments();
    unsigned int numIsoCompartments = testMCM->GetNumberOfIsotropicCompartments();
    unsigned int numNonIsoCompartments = testMCM->GetNumberOfCompartments() - numIsoCompartments;

    typedef anima::BaseCompartment::Matrix3DType TensorType;
    std::vector <ImageType::PixelType> testTensors(numNonIsoCompartments);
    std::vector <ImageType::PixelType> refTensors(numNonIsoCompartments);
    vnl_matrix <double> workLogMatrix(3,3), workLogMatrix2(3,3);
    std::vector <unsigned int> currentPermutation(numNonIsoCompartments);
    std::vector <double> testWeights(numNonIsoCompartments);
    std::vector <double> refWeights(numNonIsoCompartments);
    ImageType::PixelType zeroVec(6);
    zeroVec.Fill(0);
    zeroVec[0] = std::log(1.5e-3);
    zeroVec[2] = std::log(1.5e-3);
    zeroVec[5] = std::log(1.5e-3);

    anima::MCML2DistanceComputer::Pointer l2NormComputer = anima::MCML2DistanceComputer::New();
    // Truncates influence of stationary water, put to 0.0 to get original one
    l2NormComputer->SetLowPassGaussianSigma(2000.0);
    l2NormComputer->SetSquaredDistance(true);

    anima::LogEuclideanTensorCalculator <double>::Pointer leCalculator = anima::LogEuclideanTensorCalculator <double>::New();

    while (!maskItr.IsAtEnd())
    {
        if (maskItr.Get() == 0)
        {
            ++maskItr;
            ++l2Itr;
            ++wseItr;
            ++distItr;
            ++testItr;
            ++refItr;

            continue;
        }

        refMCM->SetModelVector(refItr.Get());
        testMCM->SetModelVector(testItr.Get());

        double isoDist = 0;
        double isoWSE = 0;
        for (unsigned int i = 0;i < numIsoCompartments;++i)
        {
            double refWeight = refMCM->GetCompartmentWeight(i);
            double testWeight = testMCM->GetCompartmentWeight(i);
            isoWSE += (refWeight - testWeight) * (refWeight - testWeight);

            if ((testWeight <= 0)||(refWeight <= 0))
                continue;

            leCalculator->GetTensorLogarithm(refMCM->GetCompartment(i)->GetDiffusionTensor().GetVnlMatrix().as_matrix(),
                                             workLogMatrix);
            leCalculator->GetTensorLogarithm(testMCM->GetCompartment(i)->GetDiffusionTensor().GetVnlMatrix().as_matrix(),
                                             workLogMatrix2);

            for (unsigned int j = 0;j < 3;++j)
                isoDist += (workLogMatrix(j,j) - workLogMatrix2(j,j)) * (workLogMatrix(j,j) - workLogMatrix2(j,j));
        }

        for (unsigned int i = numIsoCompartments;i < numCompartments;++i)
        {
            unsigned int index = i - numIsoCompartments;
            refWeights[index] = refMCM->GetCompartmentWeight(i);
            if (refWeights[index] > 0)
            {
                leCalculator->GetTensorLogarithm(refMCM->GetCompartment(i)->GetDiffusionTensor().GetVnlMatrix().as_matrix(),
                                                 workLogMatrix);
                anima::GetVectorRepresentation(workLogMatrix,refTensors[index],6,true);
            }
            else
                refTensors[index] = zeroVec;

            testWeights[index] = testMCM->GetCompartmentWeight(i);
            if (testWeights[index] > 0)
            {
                leCalculator->GetTensorLogarithm(testMCM->GetCompartment(i)->GetDiffusionTensor().GetVnlMatrix().as_matrix(),
                                                 workLogMatrix);
                anima::GetVectorRepresentation(workLogMatrix,testTensors[index],6,true);
            }
            else
                testTensors[index] = zeroVec;

            currentPermutation[index] = index;
        }

        double optimalDist = -1;
        double optimalWSE = -1;
        do
        {
            double distPermutation = 0;
            for (unsigned int i = 0;i < numNonIsoCompartments;++i)
            {
                unsigned int index = currentPermutation[i];
                double dist = 0;
                for (unsigned int j = 0;j < refTensors[i].GetSize();++j)
                    dist += (refTensors[i][j] - testTensors[index][j]) * (refTensors[i][j] - testTensors[index][j]);

                distPermutation += dist;
            }

            if ((distPermutation < optimalDist)||(optimalDist < 0))
            {
                optimalWSE = 0;
                for (unsigned int i = 0;i < numNonIsoCompartments;++i)
                    optimalWSE += (refWeights[i] - testWeights[currentPermutation[i]]) * (refWeights[i] - testWeights[currentPermutation[i]]);

                optimalDist = distPermutation;
            }
        } while(std::next_permutation(currentPermutation.begin(),currentPermutation.end()));

        double optimalWSEInd = -1;
        do
        {
            double distPermutation = 0;
            for (unsigned int i = 0;i < numNonIsoCompartments;++i)
                distPermutation += (refWeights[i] - testWeights[currentPermutation[i]]) * (refWeights[i] - testWeights[currentPermutation[i]]);

            if ((distPermutation < optimalWSEInd)||(optimalWSEInd < 0))
                optimalWSEInd = distPermutation;

        } while(std::next_permutation(currentPermutation.begin(),currentPermutation.end()));

        wseItr.Set(optimalWSE + isoWSE);
        wseIndItr.Set(optimalWSEInd + isoWSE);
        distItr.Set(optimalDist + isoDist);

        double l2Dist = l2NormComputer->ComputeDistance(refMCM,testMCM);
        l2Itr.Set(l2Dist);

        ++maskItr;
        ++l2Itr;
        ++wseItr;
        ++wseIndItr;
        ++distItr;
        ++testItr;
        ++refItr;
    }

    anima::writeImage<OutputImageType>(outWSEArg.getValue(),wseImage);
    anima::writeImage<OutputImageType>(outWSEIndArg.getValue(),wseIndImage);
    anima::writeImage<OutputImageType>(outArg.getValue(),distImage);
    anima::writeImage<OutputImageType>(outL2Arg.getValue(),l2Image);

    return EXIT_SUCCESS;
}
