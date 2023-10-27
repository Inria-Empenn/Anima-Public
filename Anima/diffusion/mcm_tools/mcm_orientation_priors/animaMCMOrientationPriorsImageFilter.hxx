#pragma once
#include "animaMCMOrientationPriorsImageFilter.h"
#include <animaSpectralClusteringFilter.h>

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

namespace anima
{

template <class TPixelType>
void
MCMOrientationPriorsImageFilter <TPixelType>
::GenerateOutputInformation()
{
    //-------------------------------
    // Creates the additional outputs
    //-------------------------------
    unsigned int prevNum = this->GetNumberOfOutputs();

    InputImageType *input = const_cast <InputImageType *> (this->GetInput(0));
    MCModelPointer tmpMCM = input->GetDescriptionModel();
    unsigned int numIsoCompartments = tmpMCM->GetNumberOfIsotropicCompartments();
    unsigned int numCompartments = tmpMCM->GetNumberOfCompartments();
    m_NumberOfAnisotropicCompartments = numCompartments - numIsoCompartments;
    unsigned int numOutputs = m_NumberOfAnisotropicCompartments + 1;
    
    this->SetNumberOfIndexedOutputs(numOutputs);
    
    for (unsigned int i = prevNum;i < numOutputs;++i)
        this->SetNthOutput(i,this->MakeOutput(i).GetPointer());
    
    //---------------------------------------------------
    // Call the superclass' implementation of this method
    //---------------------------------------------------
    Superclass::GenerateOutputInformation();

    //--------------------------------------------------
    // Set length of vectors in the output vector images
    //--------------------------------------------------
    for (unsigned int i = 0;i < m_NumberOfAnisotropicCompartments;++i)
    {
        OutputImageType *output = this->GetOutput(i);
        output->SetVectorLength(4); // AST: now fixed for Watson
    }

    OutputImageType *output = this->GetOutput(m_NumberOfAnisotropicCompartments);
    output->SetVectorLength(m_NumberOfAnisotropicCompartments);
}

template <class TPixelType>
void
MCMOrientationPriorsImageFilter <TPixelType>
::BeforeThreadedGenerateData()
{
    this->Superclass::BeforeThreadedGenerateData();

    m_ReferenceInputModels.resize(this->GetNumberOfIndexedInputs());

    InputImageType *input = const_cast <InputImageType *> (this->GetInput(0));
    MCModelPointer tmpMCM = input->GetDescriptionModel();
    unsigned int numIsoCompartments = tmpMCM->GetNumberOfIsotropicCompartments();
    unsigned int numCompartments = tmpMCM->GetNumberOfCompartments();

    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
    {
        InputImageType *input = const_cast <InputImageType *> (this->GetInput(i));
        MCModelPointer model = input->GetDescriptionModel();

        if (model->GetNumberOfIsotropicCompartments() != numIsoCompartments || model->GetNumberOfCompartments() != numCompartments)
            itkExceptionMacro("All input MCM images should have the same compartment structure.");

        m_ReferenceInputModels[i] = model;
    }
}

template <class TPixelType>
void
MCMOrientationPriorsImageFilter <TPixelType>
::DynamicThreadedGenerateData(const InputRegionType &region)
{
    unsigned int numInputs = this->GetNumberOfIndexedInputs();
    unsigned int numOutputs = this->GetNumberOfIndexedOutputs();

    using InputIteratorType = itk::ImageRegionConstIterator <InputImageType>;
    std::vector<InputIteratorType> inputIterators(numInputs);
    for (unsigned int i = 0; i < numInputs; ++i)
        inputIterators[i] = InputIteratorType(this->GetInput(i), region);

    using OutputIteratorType = itk::ImageRegionIterator <OutputImageType>;
    std::vector<OutputIteratorType> outputIterators(numOutputs);
    for (unsigned int i = 0; i < numInputs; ++i)
        outputIterators[i] = OutputIteratorType(this->GetOutput(i), region);

    using MaskIteratorType = itk::ImageRegionConstIterator <MaskImageType>;
    std::vector<MaskIteratorType> maskIterators;
    if (m_MaskImages.size() == numInputs)
    {
        for (unsigned int i = 0; i < numInputs; ++i)
            maskIterators.push_back(MaskIteratorType(m_MaskImages[i], region));
    }
    unsigned int numMasks = maskIterators.size();

    std::vector<MCModelPointer> workInputModels(numInputs);
    for (unsigned int i = 0;i < numInputs;++i)
        workInputModels[i] = m_ReferenceInputModels[i]->Clone();
    
    using ClusterFilterType = anima::SpectralClusteringFilter<double>;
    using OrientationVectorType = itk::Vector<double,3>;
    ClusterFilterType clusterFilter;
    ClusterFilterType::MatrixType squaredDistanceMatrix;
    std::vector<OrientationVectorType> workAllOrientations(numInputs * m_NumberOfAnisotropicCompartments);
    std::vector<OrientationVectorType> workInputOrientations;
    std::vector<double> workAllWeights(numInputs * m_NumberOfAnisotropicCompartments);
    vnl_matrix<double> workInputWeights;
    std::vector<bool> usefulInputsForWeights(numInputs, false);
    std::vector<unsigned int> workAllMemberships(numInputs * m_NumberOfAnisotropicCompartments);
    OrientationVectorType workSphericalOrientation, workCartesianOrientation;
    workSphericalOrientation[2] = 1.0;
    std::vector<unsigned int> clusterMembers;
    // anima::WatsonDistribution watsonDistribution;
    // anima::DirichletDistribution dirichletDistribution;

    std::vector< itk::VariableLengthVector<double> > outputOrientationValues(m_NumberOfAnisotropicCompartments);
    for (unsigned int i = 0;i < m_NumberOfAnisotropicCompartments;++i)
        outputOrientationValues[i].SetSize(4);
    itk::VariableLengthVector<double> outputWeightValues;
    outputWeightValues.SetSize(m_NumberOfAnisotropicCompartments);

    while (!inputIterators[0].IsAtEnd())
    {
        std::fill(workInputWeights.begin(),workInputWeights.end(),0.0);

        unsigned int numUsefulInputs = 0;
        for (unsigned int i = 0;i < numInputs;++i)
        {
            if (isZero(inputIterators[i].Get()))
                continue;

            if (numMasks == numInputs)
            {
                if (maskIterators[i].Get() == 0)
                    continue;
            }

            workInputModels[i]->SetModelVector(inputIterators[i].Get());

            ++numUsefulInputs;
        }

        if (numUsefulInputs == 0)
        {
            for (unsigned int i = 0;i < m_NumberOfAnisotropicCompartments;++i)
            {
                outputOrientationValues[i].Fill(0.0);
                outputIterators[i].Set(outputOrientationValues[i]);
                ++outputIterators[i];
            }

            outputWeightValues.Fill(0.0);
            outputIterators[m_NumberOfAnisotropicCompartments].Set(outputWeightValues);
            ++outputIterators[m_NumberOfAnisotropicCompartments];

            for (unsigned int i = 0;i < numInputs;++i)
                ++inputIterators[i];
            
            if (numMasks == numInputs)
            {
                for (unsigned int i = 0;i < numMasks;++i)
                    ++maskIterators[i];
            }

            continue;
        }

        // Actual work comes here

        // Grab list of all orientations, memberships and weights from all MCM models
        for (unsigned int i = 0;i < numInputs;++i)
        {
            unsigned int numIsoCompartments = workInputModels[i]->GetNumberOfIsotropicCompartments();
            
            double sumOfAnisotropicWeights = 0.0;
            for (unsigned int j = 0;j < m_NumberOfAnisotropicCompartments;++j)
            {
                workSphericalOrientation[0] = workInputModels[i]->GetCompartment(j + numIsoCompartments)->GetOrientationTheta();
                workSphericalOrientation[1] = workInputModels[i]->GetCompartment(j + numIsoCompartments)->GetOrientationPhi();
                anima::TransformSphericalToCartesianCoordinates(workSphericalOrientation, workCartesianOrientation);
                workAllOrientations[i * m_NumberOfAnisotropicCompartments + j] = workCartesianOrientation;
                workAllMemberships[i * m_NumberOfAnisotropicCompartments + j] = i;
                sumOfAnisotropicWeights += workInputModels[i]->GetCompartmentWeight(j + numIsoCompartments);
            }

            for (unsigned int j = 0;j < m_NumberOfAnisotropicCompartments;++j)
                workAllWeights[i * m_NumberOfAnisotropicCompartments + j] = workInputModels[i]->GetCompartmentWeight(j + numIsoCompartments) / sumOfAnisotropicWeights;
        }

        // Compute squared distance matrix between orientations
        unsigned int numOrientations = workAllOrientations.size();
        squaredDistanceMatrix.set_size(numOrientations, numOrientations);
        squaredDistanceMatrix.fill_diagonal(0.0);
        for (unsigned int i = 0;i < numOrientations - 1;++i)
        {
            workCartesianOrientation = workAllOrientations[i];
            for (unsigned int j = i + 1;j < numOrientations;++j)
            {
                double cosValue = anima::ComputeScalarProduct(workAllOrientations[j], workCartesianOrientation);
                cosValue = std::abs(cosValue);
                if (cosValue > 1.0)
                    cosValue = 1.0;
                double distanceValue = std::acos(cosValue) / M_PI;
                squaredDistanceMatrix.put(i, j, distanceValue * distanceValue);
                if (j != i)
                    squaredDistanceMatrix.put(i, j, distanceValue * distanceValue);
            }
        }

        // Perform clustering based on orientations
        clusterFilter.SetInputData(squaredDistanceMatrix);
        clusterFilter.SetNbClass(m_NumberOfAnisotropicCompartments);
        clusterFilter.Update();

        // Characterize each cluster by
        // - A Watson distribution for orientations,
        // - A Beta distribution for the weight.
        for (unsigned int i = 0;i < m_NumberOfAnisotropicCompartments;++i)
        {
            clusterMembers = clusterFilter.GetClassMembers(i);
            unsigned int clusterSize = clusterMembers.size();
            workInputOrientations.resize(clusterSize);
            for (unsigned int j = 0;j < clusterSize;++j)
                workInputOrientations[j] = workAllOrientations[clusterMembers[j]];
                
            // watsonDistribution.Fit(workInputOrientations, "mle");
            // workCartesianOrientation = watsonDistribution.GetMeanAxis();
            // for (unsigned int j = 0;j < 3;++j)
            //     outputOrientationValues[i][j] = workCartesianOrientation[j];
            // outputOrientationValues[i][3] = watsonDistribution.GetConcentrationParameter();
            outputIterators[i].Set(outputOrientationValues[i]);
            ++outputIterators[i];
        }

        std::fill(usefulInputsForWeights.begin(), usefulInputsForWeights.end(), false);
        unsigned int nbUsefulInputsForWeights = 0;
        for (unsigned int i = 0;i < numInputs;++i)
        {
            unsigned int numCompartments = 0;
            for (unsigned int j = 0;j < m_NumberOfAnisotropicCompartments;++j)
            {
                clusterMembers = clusterFilter.GetClassMembers(j);
                for (unsigned int k = 0;k < clusterMembers.size();++k)
                {
                    if (workAllMemberships[clusterMembers[k]] == i)
                    {
                        numCompartments++;
                        break;
                    }
                }
            }

            if (numCompartments == m_NumberOfAnisotropicCompartments)
            {
                usefulInputsForWeights[i] = true;
                nbUsefulInputsForWeights++;
            }
        }

        workInputWeights.set_size(nbUsefulInputsForWeights, m_NumberOfAnisotropicCompartments);
        unsigned int pos = 0;
        for (unsigned int i = 0;i < numInputs;++i)
        {
            unsigned int numCompartments = 0;
            for (unsigned int j = 0;j < m_NumberOfAnisotropicCompartments;++j)
            {
                clusterMembers = clusterFilter.GetClassMembers(j);
                for (unsigned int k = 0;k < clusterMembers.size();++k)
                {
                    if (workAllMemberships[clusterMembers[k]] == i)
                    {
                        workInputWeights.put(pos, j, workAllWeights[clusterMembers[k]]);
                        break;
                    }
                }
            }

            if (usefulInputsForWeights[i])
                pos++;
        }

        // dirichletDistribution.Fit(workInputWeights, "mle");
        // workAllWeights = dirichletDistribution.GetConcentrationParameters();
        // for (unsigned int i = 0;i < m_NumberOfAnisotropicCompartments;++i)
        //     outputWeightValues[i] = workAllWeights[i];
        outputIterators[m_NumberOfAnisotropicCompartments].Set(outputWeightValues);
        ++outputIterators[m_NumberOfAnisotropicCompartments];

        for (unsigned int i = 0;i < numInputs;++i)
            ++inputIterators[i];
        
        if (numMasks == numInputs)
        {
            for (unsigned int i = 0;i < numMasks;++i)
                ++maskIterators[i];
        }
    }
}

} // end namespace anima
