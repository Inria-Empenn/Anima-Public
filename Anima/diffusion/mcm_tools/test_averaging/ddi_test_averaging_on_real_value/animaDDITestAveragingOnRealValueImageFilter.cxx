#include "animaDDITestAveragingOnRealValueImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>

#include <animaDDIAveragingTools.h>
#include <animaMCMWeightedAverager.h>

namespace anima
{

void DDITestAveragingOnRealValueImageFilter::BeforeThreadedGenerateData()
{
    this->Superclass::BeforeThreadedGenerateData();

    InputImageType *input = const_cast <InputImageType *> (this->GetInput());
    MCModelPointer model = input->GetDescriptionModel();
    unsigned int numIsoCompartments = model->GetNumberOfIsotropicCompartments();

    if (model->GetNumberOfCompartments() > numIsoCompartments)
    {
        if (model->GetCompartment(numIsoCompartments)->GetCompartmentType() != anima::DDI)
            itkExceptionMacro("DDI averaging only supports models with DDI compartments");
    }

    m_ReferenceInputModel = model;

    this->GetOutput()->SetDescriptionModel(m_ReferenceOutputModel);
}


void DDITestAveragingOnRealValueImageFilter::SetReferenceOutputModel(MCModelPointer &model)
{
    if (model->GetNumberOfCompartments() > model->GetNumberOfIsotropicCompartments())
    {
        if (model->GetCompartment(model->GetNumberOfIsotropicCompartments())->GetCompartmentType() != anima::DDI)
            itkExceptionMacro("DDI averaging only supports models with DDI compartments");
    }

    m_ReferenceOutputModel = model;
}


void DDITestAveragingOnRealValueImageFilter::DynamicThreadedGenerateData(const InputRegionType &region)
{
    typedef itk::ImageRegionConstIteratorWithIndex <InputImageType> InputIteratorType;
    typedef itk::ImageRegionIterator<OutputImageType> OutputIteratorType;

    InputIteratorType inputIterator(this->GetInput(), region);
    OutputIteratorType outputIterator (this->GetOutput(), region);

    InputImageType::SizeType totalSize = region.GetSize();

    typedef anima::MCMWeightedAverager::Pointer MCMAveragerPointer;
    MCMAveragerPointer mcmAverager = anima::MCMWeightedAverager::New();
    mcmAverager->SetDDIInterpolationMethod(m_Method);
    mcmAverager->SetOutputModel(m_ReferenceOutputModel);

    PixelType voxelOutputValue(m_ReferenceOutputModel->GetSize());

    while (!inputIterator.IsAtEnd())
    {
        voxelOutputValue.Fill(0.0);
        InputIndexType inputIndex = inputIterator.GetIndex();
        std::vector<InputIndexType> variableIndex;
        std::vector<double> wPosition;
        unsigned int nbVoxels = 0;

        // Compute which voxels we interpolate and their weights.
        anima::ComputeVoxelWeights(inputIndex, variableIndex, wPosition, totalSize, m_Step, nbVoxels);

        if (nbVoxels == 0)
        {
            outputIterator.Set(inputIterator.Get());

            ++inputIterator;
            ++outputIterator;

            continue;
        }

        std::vector <double> workInputWeights(nbVoxels,0.0);

        std::vector <MCModelPointer> workInputModels(nbVoxels);
        for (unsigned int i = 0;i < nbVoxels;++i)
            workInputModels[i] = m_ReferenceInputModel->Clone();

        //Pull all fascicles (and free water) with corresponding weights in one vector
        for (int i = 0;i < nbVoxels;++i)
        {
            PixelType ddiValue = this->GetInput()->GetPixel(variableIndex[i]);
            workInputModels[i]->SetModelVector(ddiValue);
            workInputWeights[i] = 1.0;
        }

        mcmAverager->SetInputModels(workInputModels);
        mcmAverager->SetInputWeights(workInputWeights);

        // weights are normalized inside averager
        mcmAverager->Update();

        voxelOutputValue = mcmAverager->GetOutputModel()->GetModelVector();

        //Set output value then upgrade iterators
        outputIterator.Set(voxelOutputValue);

        ++inputIterator;
        ++outputIterator;
    }
}

} // end namespace anima
