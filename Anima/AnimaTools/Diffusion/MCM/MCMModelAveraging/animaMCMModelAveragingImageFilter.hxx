#pragma once
#include "animaMCMModelAveragingImageFilter.h"
#include <animaMCMImageSimplifier.h>

#include <animaMultiCompartmentModelCreator.h>
#include <animaVectorOperations.h>

#include <animaModularityClusteringFilter.h>
#include <animaHyperbolicFunctions.h>

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

namespace anima
{

template <class PixelScalarType> const double MCMModelAveragingImageFilter<PixelScalarType>::m_ZeroThreshold = 1e-6;

template <class PixelScalarType>
void
MCMModelAveragingImageFilter<PixelScalarType>
::SetAICcVolume(unsigned int i, ScalarImageType *vol)
{
    if (i == m_AICcVolumes.size())
        m_AICcVolumes.push_back(vol);
    else if (i > m_AICcVolumes.size())
    {
        itkExceptionMacro("Trying to add a non contiguous AICc volume... Add AICc volumes contiguously (0,1,2,3,...)...");
    }
    else
        m_AICcVolumes[i] = vol;
}

template <class PixelScalarType>
void
MCMModelAveragingImageFilter<PixelScalarType>
::SetB0Volume(unsigned int i, ScalarImageType *vol)
{
    if (i == m_B0Volumes.size())
        m_B0Volumes.push_back(vol);
    else if (i > m_B0Volumes.size())
    {
        itkExceptionMacro("Trying to add a non contiguous B0 volume... Add B0 volumes contiguously (0,1,2,3,...)...");
    }
    else
        m_B0Volumes[i] = vol;
}

template <class PixelScalarType>
void
MCMModelAveragingImageFilter<PixelScalarType>
::SetNoiseVolume(unsigned int i, ScalarImageType *vol)
{
    if (i == m_NoiseVolumes.size())
        m_NoiseVolumes.push_back(vol);
    else if (i > m_NoiseVolumes.size())
    {
        itkExceptionMacro("Trying to add a non contiguous noise volume... Add noise volumes contiguously (0,1,2,3,...)...");
    }
    else
        m_NoiseVolumes[i] = vol;
}

template <class PixelScalarType>
void
MCMModelAveragingImageFilter<PixelScalarType>
::GenerateOutputInformation()
{
    // Override the method in itkImageSource, so we can set the vector length of
    // the output itk::VectorImage
    
    this->Superclass::GenerateOutputInformation();
    this->InitializeReferenceOutputModel();
    
    OutputImageType *output = this->GetOutput();
    output->SetVectorLength(m_ReferenceOutputModel->GetSize());
}

template <class PixelScalarType>
void
MCMModelAveragingImageFilter<PixelScalarType>
::InitializeReferenceOutputModel()
{
    if (m_ReferenceModels.size() == 0)
    {
        m_ReferenceModels.resize(this->GetNumberOfIndexedInputs());
        for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
        {
            InputImageType *input = const_cast <InputImageType *> (this->GetInput(i));
            m_ReferenceModels[i] = input->GetDescriptionModel();
        }
    }
    
    unsigned int numberOfCombinations = 1;
    
    m_WorkNonFreeWaterCorrespondences.clear();
    m_NumberOfNonFreeWaterCompartments.clear();
    m_NumberOfIsotropicCompartments = 0;
    
    // We assume that all non free water compartments are of the same type
    anima::DiffusionModelCompartmentType compartmentType;
    bool modelWithIRW = false;
    bool modelWithSW = false;
    bool modelWithFW = false;
    bool modelWithStanisz = false;

    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
    {
        unsigned int numberOfCompartments = m_ReferenceModels[i]->GetNumberOfCompartments();
        
        if (i == 0)
        {
            m_NumberOfIsotropicCompartments = m_ReferenceModels[i]->GetNumberOfIsotropicCompartments();
            for (unsigned int j = 0;j < m_NumberOfIsotropicCompartments;++j)
            {
                anima::DiffusionModelCompartmentType isoCompType = m_ReferenceModels[i]->GetCompartment(j)->GetCompartmentType();
                
                switch (isoCompType)
                {
                    case anima::FreeWater:
                        modelWithFW = true;
                        break;
                    case anima::StationaryWater:
                        modelWithSW = true;
                        break;
                    case anima::Stanisz:
                        modelWithStanisz = true;
                        break;
                    case anima::IsotropicRestrictedWater:
                    default:
                        modelWithIRW = true;
                        break;
                }
            }
        }
        
        if (m_ReferenceModels[i]->GetNumberOfIsotropicCompartments() != m_NumberOfIsotropicCompartments)
            itkExceptionMacro("All input models should have the same number of isotropic compartments.");
        
        unsigned int numberOfNonFreeWaterCompartments = numberOfCompartments - m_NumberOfIsotropicCompartments;
        
        if (numberOfNonFreeWaterCompartments > 0)
        {
            compartmentType = m_ReferenceModels[i]->GetCompartment(m_NumberOfIsotropicCompartments)->GetCompartmentType();
            
            numberOfCombinations *= numberOfNonFreeWaterCompartments;
            m_WorkNonFreeWaterCorrespondences.push_back(i);
            m_NumberOfNonFreeWaterCompartments.push_back(numberOfNonFreeWaterCompartments);
        }
    }
    
    typedef anima::MultiCompartmentModelCreator MCMCreatorType;
    MCMCreatorType mcmCreator;
    
    mcmCreator.SetModelWithFreeWaterComponent(modelWithFW);
    mcmCreator.SetModelWithStationaryWaterComponent(modelWithSW);
    mcmCreator.SetModelWithRestrictedWaterComponent(modelWithIRW);
    mcmCreator.SetModelWithStaniszComponent(modelWithStanisz);
    
    mcmCreator.SetCompartmentType(compartmentType);
    mcmCreator.SetNumberOfCompartments(3 * numberOfCombinations);
    
    m_ReferenceOutputModel = mcmCreator.GetNewMultiCompartmentModel();
    this->GetOutput()->SetDescriptionModel(m_ReferenceOutputModel);
}

template <class PixelScalarType>
void
MCMModelAveragingImageFilter<PixelScalarType>
::BeforeThreadedGenerateData()
{
    Superclass::BeforeThreadedGenerateData();
    
    if (m_AICcVolumes.size() != this->GetNumberOfIndexedInputs())
    {
        std::string error("There should be the same number of input images and input AICc volumes... ");
        error += m_AICcVolumes.size();
        error += " ";
        error += this->GetNumberOfIndexedInputs();
        throw itk::ExceptionObject(__FILE__, __LINE__,error,ITK_LOCATION);
    }
    
    if (!m_ReferenceOutputModel)
        this->InitializeReferenceOutputModel();
    
    // Create Mose map
    m_MoseMap = MoseImageType::New();
    m_MoseMap->Initialize();
    m_MoseMap->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
    m_MoseMap->SetSpacing (this->GetInput(0)->GetSpacing());
    m_MoseMap->SetOrigin (this->GetInput(0)->GetOrigin());
    m_MoseMap->SetDirection (this->GetInput(0)->GetDirection());
    m_MoseMap->Allocate();
    m_MoseMap->FillBuffer(0);
    
    if (m_B0Volumes.size() != 0)
    {
        m_OutputB0Volume = ScalarImageType::New();
        m_OutputB0Volume->Initialize();
        m_OutputB0Volume->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
        m_OutputB0Volume->SetSpacing (this->GetInput(0)->GetSpacing());
        m_OutputB0Volume->SetOrigin (this->GetInput(0)->GetOrigin());
        m_OutputB0Volume->SetDirection (this->GetInput(0)->GetDirection());
        m_OutputB0Volume->Allocate();
        m_OutputB0Volume->FillBuffer(0);
    }
    
    if (m_NoiseVolumes.size() != 0)
    {
        m_OutputNoiseVolume = ScalarImageType::New();
        m_OutputNoiseVolume->Initialize();
        m_OutputNoiseVolume->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
        m_OutputNoiseVolume->SetSpacing (this->GetInput(0)->GetSpacing());
        m_OutputNoiseVolume->SetOrigin (this->GetInput(0)->GetOrigin());
        m_OutputNoiseVolume->SetDirection (this->GetInput(0)->GetDirection());
        m_OutputNoiseVolume->Allocate();
        m_OutputNoiseVolume->FillBuffer(0);
    }
}

template <class PixelScalarType>
void
MCMModelAveragingImageFilter<PixelScalarType>
::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    // Iterator definitions
    typedef itk::ImageRegionConstIterator <InputImageType> ImageIteratorType;
    
    unsigned int numInputs = this->GetNumberOfIndexedInputs();
    std::vector <ImageIteratorType> inIterators(numInputs);
    for (unsigned int i = 0;i < numInputs;++i)
        inIterators[i] = ImageIteratorType(this->GetInput(i),outputRegionForThread);
    
    typedef itk::ImageRegionIterator <OutputImageType> OutImageIteratorType;
    OutImageIteratorType outIterator(this->GetOutput(),outputRegionForThread);
    
    typedef itk::ImageRegionIterator <MoseImageType> MoseImageIteratorType;
    MoseImageIteratorType moseIterator(m_MoseMap,outputRegionForThread);
    
    typedef itk::ImageRegionIterator <ScalarImageType> ScalarImageIteratorType;
    std::vector <ScalarImageIteratorType> aiccIterators(numInputs);
    for (unsigned int i = 0;i < numInputs;++i)
        aiccIterators[i] = ScalarImageIteratorType(m_AICcVolumes[i],outputRegionForThread);
    
    std::vector <ScalarImageIteratorType> b0Iterators(numInputs);
    ScalarImageIteratorType outputB0Iterator;
    if (m_OutputB0Volume)
    {
        for (unsigned int i = 0;i < numInputs;++i)
            b0Iterators[i] = ScalarImageIteratorType(m_B0Volumes[i],outputRegionForThread);
        
        outputB0Iterator = ScalarImageIteratorType(m_OutputB0Volume,outputRegionForThread);
    }
    
    std::vector <ScalarImageIteratorType> noiseIterators(numInputs);
    ScalarImageIteratorType outputNoiseIterator;
    if (m_OutputNoiseVolume)
    {
        for (unsigned int i = 0;i < numInputs;++i)
            noiseIterators[i] = ScalarImageIteratorType(m_NoiseVolumes[i],outputRegionForThread);
        
        outputNoiseIterator = ScalarImageIteratorType(m_OutputNoiseVolume,outputRegionForThread);
    }
    
    MCModelPointer referenceOutputModel = m_ReferenceOutputModel->Clone();
    MCModelPointer outputModel = m_ReferenceOutputModel->Clone();
    
    std::vector <MCModelPointer> referenceModels(m_ReferenceModels.size());
    for (unsigned int i = 0;i < referenceModels.size();++i)
        referenceModels[i] = m_ReferenceModels[i]->Clone();
    
    // Initialize output vector
    unsigned int outVectorSize = referenceOutputModel->GetSize();
    OutputPixelType vecImage(outVectorSize);
    OutputPixelType resVecImage(outVectorSize);
    BaseCompartmentType::ModelOutputVectorType resVecModel(outVectorSize);
    
    // Data for eigen analysis
    MatrixType eigVecs;
    VectorType eigVals;
    
    // Actual work data
    WeightDataType modelAICs(numInputs,0);
    WeightDataType modelAICWeights(numInputs,0);
    
    unsigned int pairingVectorSize = m_WorkNonFreeWaterCorrespondences.size();
    unsigned int lastModelDimension = m_NumberOfNonFreeWaterCompartments[pairingVectorSize - 1];
    std::vector <unsigned int> modelPairingVector(m_WorkNonFreeWaterCorrespondences.size(),0);
    BaseCompartmentType::ListType outputCompartmentWeights(referenceOutputModel->GetNumberOfCompartments(),0);
    vnl_matrix <double> connectivityMatrix, distanceMatrix;
    MatrixType dcmData;
    VectorType fascicleDirection, fascicleDirection2;
    std::vector <unsigned int> clusterMembers;
    
    while (!outIterator.IsAtEnd())
    {
        resVecImage.Fill(0.0);
        resVecModel.Fill(0.0);
        
        bool uselessVoxel = true;
        for (unsigned int i = 0;i < numInputs;++i)
        {
            if (aiccIterators[i].Get() != 0)
            {
                uselessVoxel = false;
                break;
            }
        }
        
        if (uselessVoxel)
        {
            outIterator.Set(resVecImage);
            ++outIterator;
            
            moseIterator.Set(0);
            ++moseIterator;
            
            for (unsigned int i = 0;i < numInputs;++i)
            {
                ++inIterators[i];
                ++aiccIterators[i];
            }
            
            if (m_OutputB0Volume)
            {
                for (unsigned int i = 0;i < numInputs;++i)
                    ++b0Iterators[i];
                
                ++outputB0Iterator;
            }
            
            if (m_OutputNoiseVolume)
            {
                for (unsigned int i = 0;i < numInputs;++i)
                    ++noiseIterators[i];
                
                ++outputNoiseIterator;
            }

            this->IncrementNumberOfProcessedPoints();
            continue;
        }
        
        std::fill(outputCompartmentWeights.begin(),outputCompartmentWeights.end(),0);
        outputModel->SetModelVector(resVecModel);
        referenceOutputModel->SetModelVector(resVecModel);
        
        for (unsigned int i = 0;i < numInputs;++i)
        {
            vecImage = inIterators[i].Get();
            referenceModels[i]->SetModelVector(vecImage);
            modelAICs[i] = aiccIterators[i].Get();
        }
        
        modelAICWeights = this->GetAkaikeWeights(modelAICs);
        
        // Start by easy and optional things
        if (m_OutputB0Volume)
        {
            double outputB0Value = 0;
            for (unsigned int i = 0;i < numInputs;++i)
                outputB0Value += modelAICWeights[i] * b0Iterators[i].Get();
            
            outputB0Iterator.Set(outputB0Value);
        }
        
        if (m_OutputNoiseVolume)
        {
            double outputNoiseValue = 0;
            for (unsigned int i = 0;i < numInputs;++i)
                outputNoiseValue += modelAICWeights[i] * noiseIterators[i].Get();
            
            outputNoiseIterator.Set(outputNoiseValue);
        }
        
        double sumIsoWeights = 0;
        // Perform free water model averaging, filling directly output model here
        if (m_NumberOfIsotropicCompartments > 0)
        {
            for (unsigned int i = 0;i < m_NumberOfIsotropicCompartments;++i)
            {
                double axialDiffusivity = 0;
                double tissueRadius = 0;
                double averageIsoWeight = 0;
                
                double sumWeights = 0;
                for (unsigned int j = 0;j < numInputs;++j)
                {
                    double weight = modelAICWeights[j];
                    double compWeight = referenceModels[j]->GetCompartmentWeight(i);
                    if (compWeight == 0)
                        continue;
                    
                    axialDiffusivity += weight * referenceModels[j]->GetCompartment(i)->GetAxialDiffusivity();
                    tissueRadius += weight * referenceModels[j]->GetCompartment(i)->GetTissueRadius();
                    averageIsoWeight += weight * compWeight;
                    sumWeights += weight;
                }
                
                if (sumWeights > 0)
                {
                    outputModel->GetCompartment(i)->SetAxialDiffusivity(axialDiffusivity / sumWeights);
                    outputModel->GetCompartment(i)->SetTissueRadius(tissueRadius / sumWeights);
                }

                outputCompartmentWeights[i] = averageIsoWeight;
                sumIsoWeights += averageIsoWeight;
            }
        }
        
        unsigned int modelOutputPosition = m_NumberOfIsotropicCompartments;
        
        if (sumIsoWeights > 1.0 - m_WeightThreshold)
        {
            for (unsigned int i = 0;i < m_NumberOfIsotropicCompartments;++i)
                outputCompartmentWeights[i] /= sumIsoWeights;
            
            sumIsoWeights = 1.0;
        }
        else
        {
            // Now directional models averaging
            std::fill(modelPairingVector.begin(),modelPairingVector.end(),0);
            unsigned int numberOfPossibleCombinations = std::tgamma(lastModelDimension+1);
            
            while (modelPairingVector[pairingVectorSize - 1] < lastModelDimension)
            {
                // Compute combination
                dcmData.fill(0);
                double sumWeights = 0;
                double averagePerpendicularAngle = 0;
                double averageAxialDiffusivity = 0;
                double averageTissueRadius = 0;
                double averageRadialDiffusivity1 = 0;
                double averageRadialDiffusivity2 = 0;
                double averageOrientationConcentration = 0;
                double averageExtraAxonalFraction = 0;
                double averageCompartmentWeight = 0;
                
                for (unsigned int i = 0;i < modelPairingVector.size();++i)
                {
                    unsigned int modelIndex = m_WorkNonFreeWaterCorrespondences[i];
                    unsigned int compartmentIndex = m_NumberOfIsotropicCompartments + modelPairingVector[i];
                    
                    double weight = modelAICWeights[modelIndex];
                    double compWeight = referenceModels[modelIndex]->GetCompartmentWeight(compartmentIndex) * modelIndex / numberOfPossibleCombinations;
                    if ((weight == 0.0)||(compWeight == 0.0))
                        continue;
                    
                    sumWeights += weight;
                    BaseCompartmentType *modelCompartment = referenceModels[modelIndex]->GetCompartment(compartmentIndex);
                    anima::TransformSphericalToCartesianCoordinates(modelCompartment->GetOrientationTheta(),modelCompartment->GetOrientationPhi(),1.0,fascicleDirection);
                    
                    for (unsigned int j = 0;j < 3;++j)
                        for (unsigned int k = j;k < 3;++k)
                        {
                            double val = weight * fascicleDirection[j] * fascicleDirection[k];
                            dcmData(j,k) += val;
                        }
                    
                    averagePerpendicularAngle += weight * modelCompartment->GetPerpendicularAngle();
                    averageAxialDiffusivity += weight * modelCompartment->GetAxialDiffusivity();
                    averageTissueRadius += weight * modelCompartment->GetTissueRadius();
                    averageRadialDiffusivity1 += weight * modelCompartment->GetRadialDiffusivity1();
                    averageRadialDiffusivity2 += weight * modelCompartment->GetRadialDiffusivity2();
                    averageOrientationConcentration += weight * modelCompartment->GetOrientationConcentration();
                    averageExtraAxonalFraction += weight * modelCompartment->GetExtraAxonalFraction();
                    averageCompartmentWeight += weight * compWeight;
                }
                
                if (sumWeights < m_ZeroThreshold)
                {
                    // Increment pairing vector and get to next combination
                    this->IncrementModelPairingVector(modelPairingVector);
                    continue;
                }
                
                for (unsigned int j = 0;j < 3;++j)
                    for (unsigned int k = j + 1;k < 3;++k)
                        dcmData(k,j) = dcmData(j,k);
                
                dcmData /= sumWeights;
                averagePerpendicularAngle /= sumWeights;
                averageAxialDiffusivity /= sumWeights;
                averageTissueRadius /= sumWeights;
                averageRadialDiffusivity1 /= sumWeights;
                averageRadialDiffusivity2 /= sumWeights;
                averageOrientationConcentration /= sumWeights;
                averageExtraAxonalFraction /= sumWeights;
                
                EigenAnalysisType eigSystem(3);
                eigSystem.SetOrderEigenValues(true);
                eigSystem.ComputeEigenValuesAndVectors(dcmData, eigVals, eigVecs);
                
                unsigned int minPosition = modelOutputPosition;
                sumWeights = 0;
                
                for (unsigned int j = 0;j < 3;++j)
                {
                    double eigVal = eigVals[j];
                    
                    if (3.0 * eigVal < 1.0)
                        continue;
                    
                    sumWeights += eigVal;
                    
                    anima::TransformCartesianToSphericalCoordinates(eigVecs.get_row(j),fascicleDirection);
                    
                    BaseCompartmentType *outputModelCompartment = referenceOutputModel->GetCompartment(modelOutputPosition);
                    
                    outputModelCompartment->SetOrientationTheta(fascicleDirection[0]);
                    outputModelCompartment->SetOrientationPhi(fascicleDirection[1]);
                    outputModelCompartment->SetPerpendicularAngle(averagePerpendicularAngle);
                    outputModelCompartment->SetTissueRadius(averageTissueRadius);
                    outputModelCompartment->SetAxialDiffusivity(averageAxialDiffusivity);
                    outputModelCompartment->SetRadialDiffusivity1(averageRadialDiffusivity1);
                    outputModelCompartment->SetRadialDiffusivity2(averageRadialDiffusivity2);
                    outputModelCompartment->SetOrientationConcentration(averageOrientationConcentration);
                    outputModelCompartment->SetExtraAxonalFraction(averageExtraAxonalFraction);
                    outputCompartmentWeights[modelOutputPosition] = averageCompartmentWeight * eigVal;
                    
                    ++modelOutputPosition;
                }
                
                unsigned int maxPosition = modelOutputPosition;
                
                while (minPosition < maxPosition)
                {
                    outputCompartmentWeights[minPosition] /= sumWeights;
                    ++minPosition;
                }
                
                // Increment pairing vector
                this->IncrementModelPairingVector(modelPairingVector);
            }
        }
        
        double checkSum = 0;
        for (unsigned int i = 0;i < outputCompartmentWeights.size();++i)
            checkSum += outputCompartmentWeights[i];
        
        if (std::abs(checkSum - 1.0) > m_ZeroThreshold)
            itkExceptionMacro("ERROR: Compartment weights do not sum to one");
        
        referenceOutputModel->SetCompartmentWeights(outputCompartmentWeights);
        
        unsigned int realNumberOfCombinations = modelOutputPosition - m_NumberOfIsotropicCompartments;
        
        if (m_SimplifyModels && realNumberOfCombinations > 1)
        {
            // Now perform output model simplification            
            // modularity clustering wants adjacency matrix
            
            connectivityMatrix.set_size(realNumberOfCombinations,realNumberOfCombinations);
            connectivityMatrix.fill(0.0);
            
            double anisotropicWeightSum = 0;
            for (unsigned int i = m_NumberOfIsotropicCompartments;i < modelOutputPosition;++i)
                anisotropicWeightSum += outputCompartmentWeights[i];
            
            for (unsigned int i = 0;i < realNumberOfCombinations;++i)
                connectivityMatrix(i,i) = outputCompartmentWeights[i+m_NumberOfIsotropicCompartments] / anisotropicWeightSum;
            
            for (unsigned int i = 0;i < realNumberOfCombinations;++i)
            {
                BaseCompartmentType *modelCompartment1 = referenceOutputModel->GetCompartment(i+m_NumberOfIsotropicCompartments);
                anima::TransformSphericalToCartesianCoordinates(modelCompartment1->GetOrientationTheta(),modelCompartment1->GetOrientationPhi(),1.0,fascicleDirection);
                
                for (unsigned int j = i+1;j < realNumberOfCombinations;++j)
                {
                    BaseCompartmentType *modelCompartment2 = referenceOutputModel->GetCompartment(j+m_NumberOfIsotropicCompartments);
                    anima::TransformSphericalToCartesianCoordinates(modelCompartment2->GetOrientationTheta(),modelCompartment2->GetOrientationPhi(),1.0,fascicleDirection2);
                    
                    double simValue = anima::ComputeScalarProduct(fascicleDirection,fascicleDirection2);
                    
                    if (m_SquaredSimilarity)
                        simValue *= simValue;
                    else
                        simValue = std::abs(simValue);
                    
                    connectivityMatrix(i,j) = simValue;
                    connectivityMatrix(j,i) = simValue;
                }
            }
            
            anima::ModularityClusteringFilter <double> modFilter;
            modFilter.SetInputData(connectivityMatrix);
            modFilter.Update();
            
            unsigned int numberOfClusters = modFilter.GetNumberOfClusters();
            
            modelOutputPosition = m_NumberOfIsotropicCompartments;
            std::fill(outputCompartmentWeights.begin()+modelOutputPosition,outputCompartmentWeights.end(),0);
            
            for (unsigned int i = 0;i < numberOfClusters;++i)
            {
                clusterMembers = modFilter.GetReverseClassMembership(i);
                unsigned int clusterSize = clusterMembers.size();
                
                BaseCompartmentType *outputModelCompartment = outputModel->GetCompartment(modelOutputPosition);
                
                dcmData.fill(0.0);
                double sumWeights = 0;
                double averagePerpendicularAngle = 0;
                double averageTissueRadius = 0;
                double averageAxialDiffusivity = 0;
                double averageRadialDiffusivity1 = 0;
                double averageRadialDiffusivity2 = 0;
                double averageOrientationConcentration = 0;
                double averageExtraAxonalFraction = 0;
                
                for (unsigned int j = 0;j < clusterSize;++j)
                {
                    unsigned int compartmentIndex = m_NumberOfIsotropicCompartments + clusterMembers[j];
                    
                    BaseCompartmentType *modelCompartment = referenceOutputModel->GetCompartment(compartmentIndex);
                    double compartmentWeight = referenceOutputModel->GetCompartmentWeight(compartmentIndex);
                    anima::TransformSphericalToCartesianCoordinates(modelCompartment->GetOrientationTheta(),modelCompartment->GetOrientationPhi(),1.0,fascicleDirection);
                    
                    // orientation averaging
                    for (unsigned int k = 0;k < 3;++k)
                    {
                        for (unsigned int l = k;l < 3;++l)
                        {
                            double t = fascicleDirection[k] * fascicleDirection[l];
                            
                            dcmData(k,l) += t;
                            
                            if (k != l)
                                dcmData(l,k) += t;
                        }
                    }
                    
                    averagePerpendicularAngle += modelCompartment->GetPerpendicularAngle();
                    averageTissueRadius += modelCompartment->GetTissueRadius();
                    averageAxialDiffusivity += modelCompartment->GetAxialDiffusivity();
                    averageRadialDiffusivity1 += modelCompartment->GetRadialDiffusivity1();
                    averageRadialDiffusivity2 += modelCompartment->GetRadialDiffusivity2();
                    averageOrientationConcentration += modelCompartment->GetOrientationConcentration();
                    averageExtraAxonalFraction += modelCompartment->GetExtraAxonalFraction();
                    sumWeights += compartmentWeight;
                }
                
                EigenAnalysisType eigSystem(3);
                eigSystem.SetOrderEigenValues(true);
                eigSystem.ComputeEigenValuesAndVectors(dcmData, eigVals, eigVecs);
                
                anima::TransformCartesianToSphericalCoordinates(eigVecs.get_row(2),fascicleDirection);
                outputModelCompartment->SetOrientationTheta(fascicleDirection[0]);
                outputModelCompartment->SetOrientationPhi(fascicleDirection[1]);
                
                averagePerpendicularAngle /= clusterSize;
                averageTissueRadius /= clusterSize;
                averageAxialDiffusivity /= clusterSize;
                averageRadialDiffusivity1 /= clusterSize;
                averageRadialDiffusivity2 /= clusterSize;
                averageOrientationConcentration /= clusterSize;
                averageExtraAxonalFraction /= clusterSize;
                
                outputModelCompartment->SetPerpendicularAngle(averagePerpendicularAngle);
                outputModelCompartment->SetTissueRadius(averageTissueRadius);
                outputModelCompartment->SetAxialDiffusivity(averageAxialDiffusivity);
                outputModelCompartment->SetRadialDiffusivity1(averageRadialDiffusivity1);
                outputModelCompartment->SetRadialDiffusivity2(averageRadialDiffusivity2);
                outputModelCompartment->SetOrientationConcentration(averageOrientationConcentration);
                outputModelCompartment->SetExtraAxonalFraction(averageExtraAxonalFraction);
                
                outputCompartmentWeights[modelOutputPosition] = sumWeights;
                
                ++modelOutputPosition;
            }
        } // end model simplification
        
        checkSum = 0;
        for (unsigned int i = 0;i < outputCompartmentWeights.size();++i)
            checkSum += outputCompartmentWeights[i];
        
        if (std::abs(checkSum - 1.0) > m_ZeroThreshold)
            itkExceptionMacro("ERROR: Compartment weights do not sum to one after simplification");
        
        outputModel->SetCompartmentWeights(outputCompartmentWeights);
        
        resVecModel = outputModel->GetModelVector();
        for (unsigned int i = 0;i < resVecModel.GetSize();++i)
            resVecImage[i] = resVecModel[i];
        
        outIterator.Set(resVecImage);
        
        moseIterator.Set(modelOutputPosition - m_NumberOfIsotropicCompartments);
        
        for (unsigned int i = 0;i < numInputs;++i)
        {
            ++inIterators[i];
            ++aiccIterators[i];
        }
        
        if (m_OutputB0Volume)
        {
            for (unsigned int i = 0;i < numInputs;++i)
                ++b0Iterators[i];
            
            ++outputB0Iterator;
        }
        
        if (m_OutputNoiseVolume)
        {
            for (unsigned int i = 0;i < numInputs;++i)
                ++noiseIterators[i];
            
            ++outputNoiseIterator;
        }
        
        ++outIterator;
        ++moseIterator;
        this->IncrementNumberOfProcessedPoints();
    }
}

template <class PixelScalarType>
void
MCMModelAveragingImageFilter<PixelScalarType>
::AfterThreadedGenerateData()
{
    unsigned int numOutputCompartments = m_ReferenceOutputModel->GetNumberOfCompartments() - m_ReferenceOutputModel->GetNumberOfIsotropicCompartments();
    if ((!m_SimplifyModels)||(numOutputCompartments <= 1))
        return;
    
    std::cout << "\nAveraging performed, now simplifying models..." << std::endl;
    typedef anima::MCMImageSimplifier <PixelScalarType> ImageSimplifierType;
    typename ImageSimplifierType::Pointer simplifier = ImageSimplifierType::New();
    
    simplifier->SetMoseVolume(m_MoseMap);
    simplifier->SetInput(this->GetOutput());
    simplifier->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    
    simplifier->Update();
    
    OutputImagePointer simplifiedImage = simplifier->GetOutput();
    simplifiedImage->DisconnectPipeline();
    
    this->SetNthOutput(0,simplifiedImage);
}

template <class PixelScalarType>
void
MCMModelAveragingImageFilter<PixelScalarType>
::IncrementModelPairingVector(std::vector <unsigned int> &modelPairingVector)
{
    unsigned int pos = 0;
    
    modelPairingVector[pos]++;
    
    while ((modelPairingVector[pos] >= m_NumberOfNonFreeWaterCompartments[pos])&&(pos < modelPairingVector.size() - 1))
    {
        modelPairingVector[pos] = 0;
        ++pos;
        modelPairingVector[pos]++;
    }
}

template <class PixelScalarType>
typename MCMModelAveragingImageFilter<PixelScalarType>::WeightDataType
MCMModelAveragingImageFilter<PixelScalarType>
::GetAkaikeWeights(const WeightDataType &aicData)
{
    unsigned int nbInputs = aicData.size();
    
    WeightDataType modelWeights(nbInputs, 0);
    
    double minAIC = aicData[0];
    for (unsigned int i = 1;i < nbInputs;++i)
    {
        double tmpAIC = aicData[i];
        
        if (tmpAIC < minAIC)
            minAIC = tmpAIC;
    }
    
    double sumModelWeights = 0;
    for (unsigned int i = 0;i < nbInputs;++i)
    {
        double tmpVal = exp((minAIC - aicData[i]) / 2.0);
        modelWeights[i] = tmpVal;
        sumModelWeights += tmpVal;
    }
    
    for (unsigned int i = 0;i < nbInputs;++i)
        modelWeights[i] /= sumModelWeights;
    
    return modelWeights;
}

} // end namespace anima
