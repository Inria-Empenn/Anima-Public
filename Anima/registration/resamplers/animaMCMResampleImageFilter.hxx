#pragma once
#include <animaMCMResampleImageFilter.h>

#include <animaMCMImage.h>
#include <animaMCMLinearInterpolateImageFunction.h>

namespace anima
{

template <typename TImageType, typename TInterpolatorPrecisionType>
unsigned int
MCMResampleImageFilter<TImageType,TInterpolatorPrecisionType>
::GetOutputVectorLength()
{
    return m_ReferenceOutputModel->GetSize();
}

/**
  * Returns a class copy with important specific variables set, not the greatest implementation especially in terms
  * of interpolators and other internal pointers copy (not done here). Use with extreme caution
  */
template <typename TImageType, typename TInterpolatorPrecisionType>
itk::LightObject::Pointer
MCMResampleImageFilter<TImageType,TInterpolatorPrecisionType>
::InternalClone() const
{
    itk::LightObject::Pointer outputPointer = Superclass::InternalClone();

    Self *castPointer = dynamic_cast <Self *> (outputPointer.GetPointer());

    if (m_ReferenceOutputModel)
        castPointer->SetReferenceOutputModel(m_ReferenceOutputModel);

    return outputPointer;
}

template <typename TImageType, typename TInterpolatorPrecisionType>
void
MCMResampleImageFilter<TImageType,TInterpolatorPrecisionType>
::BeforeThreadedGenerateData()
{
    m_WorkModels.resize(this->GetNumberOfThreads());
    for (unsigned int i = 0;i < this->GetNumberOfThreads();++i)
        m_WorkModels[i] = m_ReferenceOutputModel->Clone();

    InputImageType *output = dynamic_cast <InputImageType *> (this->GetOutput());
    output->SetDescriptionModel(m_ReferenceOutputModel);
    output->SetVectorLength(m_ReferenceOutputModel->GetSize());

    this->Superclass::BeforeThreadedGenerateData();
}

template <typename TImageType, typename TInterpolatorPrecisionType>
void
MCMResampleImageFilter<TImageType,TInterpolatorPrecisionType>
::SetReferenceOutputModel(const MCModelPointer &model)
{
    m_ReferenceOutputModel = model->Clone();
}

template <typename TImageType, typename TInterpolatorPrecisionType>
void
MCMResampleImageFilter<TImageType,TInterpolatorPrecisionType>
::InitializeInterpolator()
{
    typedef anima::MCMLinearInterpolateImageFunction <InputImageType,TInterpolatorPrecisionType> InterpolatorType;

    typename InterpolatorType::Pointer tmpInterpolator = InterpolatorType::New();
    tmpInterpolator->SetReferenceOutputModel(m_ReferenceOutputModel);
    this->SetInterpolator(tmpInterpolator.GetPointer());
}

template <typename TImageType, typename TInterpolatorPrecisionType>
void
MCMResampleImageFilter<TImageType,TInterpolatorPrecisionType>
::ReorientInterpolatedModel(const InputPixelType &interpolatedModel, vnl_matrix <double> &modelOrientationMatrix,
                            InputPixelType &orientedModel, itk::ThreadIdType threadId)
{
    m_WorkModels[threadId]->SetModelVector(interpolatedModel);

    m_WorkModels[threadId]->Reorient(modelOrientationMatrix,(!this->GetFiniteStrainReorientation()));
    orientedModel = m_WorkModels[threadId]->GetModelVector();
}

} // end namespace anima
