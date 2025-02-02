#pragma once
#include "animaMCMLinearInterpolateImageFunction.h"

#include <itkMultiThreaderBase.h>

namespace anima
{
/**
     * Define the number of neighbors
     */
template<class TInputImage, class TCoordRep>
const unsigned long
MCMLinearInterpolateImageFunction< TInputImage, TCoordRep >
::m_Neighbors = 1 << TInputImage::ImageDimension;


/**
     * Constructor
     */
template<class TInputImage, class TCoordRep>
MCMLinearInterpolateImageFunction< TInputImage, TCoordRep >
::MCMLinearInterpolateImageFunction()
{
}

template<class TInputImage, class TCoordRep>
void
MCMLinearInterpolateImageFunction< TInputImage, TCoordRep >
::SetInputImage(const InputImageType *ptr)
{
    this->Superclass::SetInputImage(ptr);

    InputImageType *input = const_cast <InputImageType *> (ptr);
    MCModelPointer model = input->GetDescriptionModel();

    if (m_MCMAveragers.size() != 0)
    {
        if (m_MCMAveragers[0].IsNotNull())
            this->TestModelsAdequation(model,m_MCMAveragers[0]->GetUntouchedOutputModel());
    }

    unsigned int numThreads = itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads();
    m_ReferenceInputModels.resize(numThreads);
    m_ReferenceInputWeights.resize(numThreads);
    for (unsigned int i = 0;i < numThreads;++i)
    {
        m_ReferenceInputModels[i].resize(m_Neighbors);
        m_ReferenceInputWeights[i].resize(m_Neighbors);
        std::fill(m_ReferenceInputWeights[i].begin(),m_ReferenceInputWeights[i].end(),0.0);
        for (unsigned int j = 0;j < m_Neighbors;++j)
            m_ReferenceInputModels[i][j] = model->Clone();
    }
}

template<class TInputImage, class TCoordRep>
void
MCMLinearInterpolateImageFunction< TInputImage, TCoordRep >
::SetReferenceOutputModel(MCModelPointer &model)
{
    if (this->GetInputImage() != 0)
        this->TestModelsAdequation(m_ReferenceInputModels[0][0],model);

    this->ResetAveragePointers(model);
}

template<class TInputImage, class TCoordRep>
void
MCMLinearInterpolateImageFunction< TInputImage, TCoordRep >
::ResetAveragePointers(MCModelPointer &model)
{
    unsigned int numThreads = itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads();
    m_MCMAveragers.resize(numThreads);

    for (unsigned int i = 0;i < numThreads;++i)
    {
        m_MCMAveragers[i] = AveragerType::New();
        m_MCMAveragers[i]->SetOutputModel(model);
    }
}

template<class TInputImage, class TCoordRep>
bool
MCMLinearInterpolateImageFunction< TInputImage, TCoordRep >
::CheckModelCompatibility(MCModelPointer &model)
{	
	unsigned int numIsoCompartments = model->GetNumberOfIsotropicCompartments();
    return  !((model->GetCompartment(numIsoCompartments)->GetCompartmentType() != anima::DDI) && (!model->GetCompartment(numIsoCompartments)->GetTensorCompatible()));
}

template<class TInputImage, class TCoordRep>
void
MCMLinearInterpolateImageFunction< TInputImage, TCoordRep >
::TestModelsAdequation(MCModelPointer &inputModel, MCModelPointer &outputModel)
{
    unsigned int numIsoCompartmentsInput = inputModel->GetNumberOfIsotropicCompartments();
    unsigned int numIsoCompartmentsOutput = outputModel->GetNumberOfIsotropicCompartments();

    if (numIsoCompartmentsInput != numIsoCompartmentsOutput)
        itkExceptionMacro("Cannot handle interpolation to model with different number of isotropic compartments.");

    for (unsigned int i = 0;i < numIsoCompartmentsInput;++i)
    {
        if (outputModel->GetCompartment(i)->GetCompartmentType() != inputModel->GetCompartment(i)->GetCompartmentType())
            itkExceptionMacro("Interpolation only of models with the same isotropic compartment types.");
    }

    if (outputModel->GetNumberOfCompartments() - numIsoCompartmentsOutput > 0)
    {
        if (outputModel->GetCompartment(numIsoCompartmentsOutput)->GetCompartmentType()
                != inputModel->GetCompartment(numIsoCompartmentsOutput)->GetCompartmentType())
            itkExceptionMacro("Interpolation of non tensor compatible models only to the same compartment type.");

        if (!this->CheckModelCompatibility(outputModel))
            itkExceptionMacro("Model not yet supported by MCM linear interpolation");
    }
}


template<class TInputImage, class TCoordRep>
void
MCMLinearInterpolateImageFunction< TInputImage, TCoordRep >
::SetSpecificAveragerParameters(unsigned int threadIndex) const
{
    typedef anima::MCMWeightedAverager InternalAveragerType;
    InternalAveragerType *castAverager = dynamic_cast <InternalAveragerType *> (this->GetAveragers()[threadIndex].GetPointer());

    castAverager->SetDDIInterpolationMethod(m_DDIInterpolationMethod);
}

template<class TInputImage, class TCoordRep>
unsigned int
MCMLinearInterpolateImageFunction< TInputImage, TCoordRep >
::GetFreeWorkIndex() const
{
    m_LockUsedModels.lock();

    unsigned int workIndex = 0;
    bool workIndexOk = false;

    while ((!workIndexOk)&&(workIndex < m_ReferenceInputModels.size()))
    {
        workIndexOk = true;
        for (unsigned int i = 0;i < m_UsedModels.size();++i)
        {
            if (m_UsedModels[i] == workIndex)
            {
                workIndexOk = false;
                break;
            }
        }

        if (!workIndexOk)
            ++workIndex;
    }

    m_UsedModels.push_back(workIndex);
    m_LockUsedModels.unlock();

    return workIndex;
}

template<class TInputImage, class TCoordRep>
void
MCMLinearInterpolateImageFunction< TInputImage, TCoordRep >
::UnlockWorkIndex(unsigned int index) const
{
    m_LockUsedModels.lock();

    for (unsigned int i = 0;i < m_UsedModels.size();++i)
    {
        if (m_UsedModels[i] == index)
        {
            m_UsedModels.erase(m_UsedModels.begin() + i);
            break;
        }
    }

    m_LockUsedModels.unlock();
}

/**
     * Evaluate at image index position
     */
template<class TInputImage, class TCoordRep>
typename MCMLinearInterpolateImageFunction< TInputImage, TCoordRep >
::OutputType
MCMLinearInterpolateImageFunction< TInputImage, TCoordRep >
::EvaluateAtContinuousIndex(const ContinuousIndexType& index) const
{
    IndexType baseIndex;
    double distance[ImageDimension], oppDistance[ImageDimension];

    for (unsigned int dim = 0; dim < ImageDimension; ++dim)
    {
        baseIndex[dim] = itk::Math::Floor< IndexValueType >( index[dim] );
        distance[dim] = index[dim] - static_cast< double >( baseIndex[dim] );
        oppDistance[dim] = 1.0 - distance[dim];
    }

    unsigned int threadIndex = this->GetFreeWorkIndex();

    OutputType voxelOutputValue(m_MCMAveragers[threadIndex]->GetOutputModelSize());
    voxelOutputValue.Fill(0);
    std::fill(m_ReferenceInputWeights[threadIndex].begin(),m_ReferenceInputWeights[threadIndex].end(),0.0);

    this->SetSpecificAveragerParameters(threadIndex);
    m_MCMAveragers[threadIndex]->ResetNumberOfOutputDirectionalCompartments();

    double totalOverlap = 0;

    unsigned int posMCM = 0;
    unsigned int maxNumCompartments = 0;

    unsigned int numIsoCompartments = m_ReferenceInputModels[threadIndex][0]->GetNumberOfIsotropicCompartments();
    unsigned int numberOfTotalInputCompartments = m_ReferenceInputModels[threadIndex][0]->GetNumberOfCompartments();

    for (unsigned int counterInput = 0;counterInput < m_Neighbors;++counterInput)
    {
        double overlap = 1.0; // fraction overlap
        unsigned int upper = counterInput; // each bit indicates upper/lower neighbour
        IndexType neighIndex;

        // get neighbor index and overlap fraction
        bool okValue = true;
        for (unsigned int dim = 0; dim < ImageDimension; ++dim)
        {
            if (upper & 1)
            {
                neighIndex[dim] = baseIndex[dim] + 1;

                if (neighIndex[dim] > this->m_EndIndex[dim])
                {
                    okValue = false;
                    break;
                }

                overlap *= distance[dim];
            }
            else
            {
                neighIndex[dim] = baseIndex[dim];

                if (neighIndex[dim] < this->m_StartIndex[dim])
                {
                    okValue = false;
                    break;
                }

                overlap *= oppDistance[dim];
            }

            upper >>= 1;
        }

        // get neighbor value only if overlap is not zero
        if ((overlap > 0) && okValue)
        {
            const PixelType input = this->GetInputImage()->GetPixel( neighIndex );
            if (isZero(input))
                continue;

            m_ReferenceInputModels[threadIndex][posMCM]->SetModelVector(input);
            m_ReferenceInputWeights[threadIndex][posMCM] = overlap;

            unsigned int numEffectiveAnisotropicCompartments = 0;
            for (unsigned int i = numIsoCompartments;i < numberOfTotalInputCompartments;++i)
            {
                double fascWeight = m_ReferenceInputModels[threadIndex][posMCM]->GetCompartmentWeight(i);
                if (fascWeight <= 0)
                    continue;

                ++numEffectiveAnisotropicCompartments;
            }

            if (numEffectiveAnisotropicCompartments > maxNumCompartments)
                maxNumCompartments = numEffectiveAnisotropicCompartments;

            totalOverlap += overlap;
            ++posMCM;
        }
    }

    if (totalOverlap < 0.5)
    {
        this->UnlockWorkIndex(threadIndex);
        voxelOutputValue.Fill(0.0);
        return voxelOutputValue;
    }

    m_MCMAveragers[threadIndex]->SetNumberOfOutputDirectionalCompartments(maxNumCompartments);
    m_MCMAveragers[threadIndex]->SetInputModels(m_ReferenceInputModels[threadIndex]);
    m_MCMAveragers[threadIndex]->SetInputWeights(m_ReferenceInputWeights[threadIndex]);

    m_MCMAveragers[threadIndex]->Update();

    voxelOutputValue = m_MCMAveragers[threadIndex]->GetOutputModel()->GetModelVector();

    this->UnlockWorkIndex(threadIndex);

    // Be sure there is no nan value in a voxel
    for (int i = 0; i < voxelOutputValue.GetSize(); ++i)
    {
        if (std::isnan(voxelOutputValue[i]))
            itkExceptionMacro("Problem, OutputVoxel has a nan value.");
    }

    return voxelOutputValue;
}

} // end namespace anima
