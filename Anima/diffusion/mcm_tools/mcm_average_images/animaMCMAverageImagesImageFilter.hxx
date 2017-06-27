#pragma once
#include "animaMCMAverageImagesImageFilter.h"

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

namespace anima
{

template <class TPixelType>
void
MCMAverageImagesImageFilter <TPixelType>
::BeforeThreadedGenerateData()
{
    this->Superclass::BeforeThreadedGenerateData();

    m_ReferenceInputModels.resize(this->GetNumberOfIndexedInputs());

    InputImageType *input = const_cast <InputImageType *> (this->GetInput(0));
    MCModelPointer tmpMCM = input->GetDescriptionModel();
    unsigned int numIsoCompartments = tmpMCM->GetNumberOfIsotropicCompartments();

    std::vector <anima::DiffusionModelCompartmentType> refIsoVec(numIsoCompartments);
    for (unsigned int j = 0;j < numIsoCompartments;++j)
        refIsoVec[j] = m_ReferenceOutputModel->GetCompartment(j)->GetCompartmentType();

    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
    {
        InputImageType *input = const_cast <InputImageType *> (this->GetInput(i));
        MCModelPointer model = input->GetDescriptionModel();

        if (model->GetNumberOfIsotropicCompartments() != numIsoCompartments)
            itkExceptionMacro("All input images should have the same isotropic compartment structure as the output model");

        for (unsigned int j = 0;j < numIsoCompartments;++j)
        {
            if (refIsoVec[j] != model->GetCompartment(j)->GetCompartmentType())
                itkExceptionMacro("All input images should have the same isotropic compartment structure as the output model");
        }

        m_ReferenceInputModels[i] = model;
    }

    this->GetOutput()->SetDescriptionModel(m_ReferenceOutputModel);
}

template <class TPixelType>
void
MCMAverageImagesImageFilter <TPixelType>
::SetReferenceOutputModel(MCModelPointer &model)
{
    m_ReferenceOutputModel = model;
}

template <class TPixelType>
typename MCMAverageImagesImageFilter <TPixelType>::MCMAveragerPointer
MCMAverageImagesImageFilter <TPixelType>
::CreateAverager()
{
    MCMAveragerPointer mcmAverager = anima::MCMWeightedAverager::New();
    mcmAverager->SetOutputModel(m_ReferenceOutputModel);

    return mcmAverager;
}

template <class TPixelType>
void
MCMAverageImagesImageFilter <TPixelType>
::ThreadedGenerateData(const InputRegionType &region, itk::ThreadIdType threadId)
{
    unsigned int numInputs = this->GetNumberOfIndexedInputs();

    typedef itk::ImageRegionConstIterator <InputImageType> InputIteratorType;
    typedef itk::ImageRegionIterator <OutputImageType> OutputIteratorType;

    std::vector<InputIteratorType> inputIterators (numInputs);

    for (unsigned int i = 0; i < numInputs; ++i)
        inputIterators[i] = InputIteratorType (this->GetInput(i), region);

    OutputIteratorType outputIterator (this->GetOutput(), region);

    MCMAveragerPointer mcmAverager = this->CreateAverager();

    std::vector <MCModelPointer> workInputModels(numInputs);
    for (unsigned int i = 0;i < numInputs;++i)
        workInputModels[i] = m_ReferenceInputModels[i]->Clone();

    std::vector <double> workInputWeights(numInputs,0.0);
    PixelType voxelOutputValue(m_ReferenceOutputModel->GetSize());

    while (!inputIterators[0].IsAtEnd())
    {
        voxelOutputValue.Fill(0.0);
        std::fill(workInputWeights.begin(),workInputWeights.end(),0.0);

        unsigned int trueInputNumber = 0;
        for (unsigned int i = 0;i < numInputs;++i)
        {
            if (isZero(inputIterators[i].Get()))
                continue;

            workInputModels[i]->SetModelVector(inputIterators[i].Get());
            workInputWeights[i] = 1.0;

            ++trueInputNumber;
        }

        if (trueInputNumber == 0)
        {
            outputIterator.Set(voxelOutputValue);
            for (unsigned int i = 0;i < numInputs;++i)
                ++inputIterators[i];
            ++outputIterator;
            continue;
        }

        mcmAverager->SetInputModels(workInputModels);
        mcmAverager->SetInputWeights(workInputWeights);

        // weights are normalized inside averager
        mcmAverager->Update();

        voxelOutputValue = mcmAverager->GetOutputModel()->GetModelVector();

        //Set output value then upgrade iterators
        outputIterator.Set(voxelOutputValue);
        for (unsigned int i = 0;i < numInputs;++i)
            ++inputIterators[i];
        ++outputIterator;
    }
}

} // end namespace anima
