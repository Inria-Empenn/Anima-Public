#pragma once
#include "animaMCMAverageImagesImageFilter.h"

#include <animaMCMWeightedAverager.h>

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
	MCMAveragerPointer outAverager = nullptr;
	
	using InternalAveragerPointer = anima::MCMDDIWeightedAverager;
	
	bool hasDDICompartment = false;
	for(int i=0; i<m_ReferenceOutputModel->GetNumberOfCompartments() && !hasDDICompartment; ++i)
	{
		hasDDICompartment = m_ReferenceOutputModel->GetCompartment(i)->GetCompartmentType() == DiffusionModelCompartmentType::DDI;
	}
	
	if(hasDDICompartment)
	{
		InternalAveragerPointer mcmAverager = InternalAveragerType::New();
		mcmAverager->SetOutputModel(this->GetReferenceOutputModel());
		mcmAverager->SetDDIInterpolationMethod(m_DDIAveragingMethod);

		outAverager = mcmAverager.GetPointer();
	}
	else
	{
		MCMAveragerPointer mcmAverager = anima::MCMWeightedAverager::New();
		mcmAverager->SetOutputModel(m_ReferenceOutputModel);
		
		outAverager = mcmAverager;
	}
	
    return outAverager;
}

template <class TPixelType>
void
MCMAverageImagesImageFilter <TPixelType>
::DynamicThreadedGenerateData(const InputRegionType &region)
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

    using MaskIteratorType = itk::ImageRegionConstIterator <MaskImageType>;
    std::vector <MaskIteratorType> maskIterators;
    if (m_MaskImages.size() == numInputs)
    {
        for (unsigned int i = 0; i < numInputs; ++i)
            maskIterators.push_back(MaskIteratorType (m_MaskImages[i], region));
    }

    unsigned int numMasks = maskIterators.size();
    while (!inputIterators[0].IsAtEnd())
    {
        voxelOutputValue.Fill(0.0);
        std::fill(workInputWeights.begin(),workInputWeights.end(),0.0);

        unsigned int trueInputNumber = 0;
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
            workInputWeights[i] = 1.0;

            ++trueInputNumber;
        }

        if (trueInputNumber == 0)
        {
            outputIterator.Set(voxelOutputValue);
            for (unsigned int i = 0;i < numInputs;++i)
                ++inputIterators[i];
            ++outputIterator;

            for (unsigned int i = 0;i < numMasks;++i)
                ++maskIterators[i];

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

        for (unsigned int i = 0;i < numMasks;++i)
            ++maskIterators[i];
    }
}

} // end namespace anima
