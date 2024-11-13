#pragma once
#include "animaMCMImageSimplifier.h"

#include <animaMultiCompartmentModelCreator.h>
#include <itkImageRegionIterator.h>

namespace anima
{

template <class PixelScalarType>
void
MCMImageSimplifier<PixelScalarType>
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
MCMImageSimplifier<PixelScalarType>
::InitializeReferenceOutputModel()
{
    // We assume that all non free water compartments are of the same type
    bool modelWithIRW = false;
    bool modelWithSW = false;
    bool modelWithFW = false;
    bool modelWithStanisz = false;

    InputImageType *inImage = const_cast <InputImageType *> (this->GetInput());
    unsigned int numIsotropicCompartments = inImage->GetDescriptionModel()->GetNumberOfIsotropicCompartments();
    for (unsigned int j = 0;j < numIsotropicCompartments;++j)
    {
        anima::DiffusionModelCompartmentType compIsoType = inImage->GetDescriptionModel()->GetCompartment(j)->GetCompartmentType();

        switch (compIsoType)
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

    anima::DiffusionModelCompartmentType compartmentType;
    if (inImage->GetDescriptionModel()->GetNumberOfCompartments() > numIsotropicCompartments)
       compartmentType = inImage->GetDescriptionModel()->GetCompartment(numIsotropicCompartments)->GetCompartmentType();
    else
        itkExceptionMacro("Only isotropic compartments, nothing to do here.")

    unsigned int maxCompartments = 0;
    itk::ImageRegionConstIterator <MoseImageType> moseItr(m_MoseMap,m_MoseMap->GetLargestPossibleRegion());

    while (!moseItr.IsAtEnd())
    {
        unsigned int moseVal = moseItr.Get();
        if (moseVal > maxCompartments)
            maxCompartments = moseVal;

        ++moseItr;
    }

    typedef anima::MultiCompartmentModelCreator MCMCreatorType;
    MCMCreatorType mcmCreator;
    mcmCreator.SetModelWithFreeWaterComponent(modelWithFW);
    mcmCreator.SetModelWithRestrictedWaterComponent(modelWithIRW);
    mcmCreator.SetModelWithStationaryWaterComponent(modelWithSW);
    mcmCreator.SetModelWithStaniszComponent(modelWithStanisz);
    mcmCreator.SetCompartmentType(compartmentType);
    mcmCreator.SetNumberOfCompartments(maxCompartments);

    m_ReferenceOutputModel = mcmCreator.GetNewMultiCompartmentModel();
    this->GetOutput()->SetDescriptionModel(m_ReferenceOutputModel);
}

template <class PixelScalarType>
void
MCMImageSimplifier<PixelScalarType>
::BeforeThreadedGenerateData()
{
    Superclass::BeforeThreadedGenerateData();

    if (!m_MoseMap)
        itkExceptionMacro("Mose map required to simplify");

    if (!m_ReferenceOutputModel)
        this->InitializeReferenceOutputModel();
}

template <class PixelScalarType>
void
MCMImageSimplifier<PixelScalarType>
::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    // Iterator definitions
    typedef itk::ImageRegionIterator <InputImageType> ImageIteratorType;

    InputImageType *inImage = const_cast <InputImageType *> (this->GetInput());
    ImageIteratorType inItr(inImage,this->GetInput()->GetLargestPossibleRegion());
    ImageIteratorType outItr(this->GetOutput(),this->GetInput()->GetLargestPossibleRegion());

    MCModelPointer inputModel = inImage->GetDescriptionModel()->Clone();
    MCModelPointer outputModel = m_ReferenceOutputModel->Clone();

    OutputPixelType inputVector, outputVector;

    MCModelType::ListType weights(outputModel->GetNumberOfCompartments());
    unsigned int numIsoCompartments = outputModel->GetNumberOfIsotropicCompartments();

    while (!outItr.IsAtEnd())
    {
        inputVector = inItr.Get();

        inputModel->SetModelVector(inputVector);
        std::fill(weights.begin(),weights.end(),0.0);

        for (unsigned int i = 0;i < numIsoCompartments;++i)
        {
            weights[i] = inputModel->GetCompartmentWeight(i);
            outputModel->GetCompartment(i)->SetCompartmentVector(inputModel->GetCompartment(i)->GetCompartmentVector());
        }

        unsigned int pos = numIsoCompartments;
        for (unsigned int i = numIsoCompartments;i < outputModel->GetNumberOfCompartments();++i)
        {
            double compWeight = inputModel->GetCompartmentWeight(i);
            if (compWeight <= 0)
                continue;

            weights[pos] = compWeight;
            outputModel->GetCompartment(pos)->SetCompartmentVector(inputModel->GetCompartment(i)->GetCompartmentVector());
            ++pos;
        }

        outputModel->SetCompartmentWeights(weights);
        outputVector = outputModel->GetModelVector();

        outItr.Set(outputVector);

        ++outItr;
        ++inItr;
        this->IncrementNumberOfProcessedPoints();
    }
}

} // end namespace anima
