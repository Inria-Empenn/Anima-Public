#pragma once
#include <animaFiniteStrainTensorResampleImageFilter.h>

#include <animaBaseTensorTools.h>

namespace anima
{

template <typename TImageType, typename TInterpolatorPrecisionType>
void
FiniteStrainTensorResampleImageFilter<TImageType, TInterpolatorPrecisionType>
::BeforeThreadedGenerateData()
{
    this->Superclass::BeforeThreadedGenerateData();

    m_WorkMats.resize(this->GetNumberOfThreads());
    m_TmpTensors.resize(this->GetNumberOfThreads());

    for (unsigned int i = 0;i < this->GetNumberOfThreads();++i)
    {
        m_WorkMats[i].set_size(m_TensorDimension,m_TensorDimension);
        m_TmpTensors[i].set_size(m_TensorDimension,m_TensorDimension);
    }
}

template <typename TImageType, typename TInterpolatorPrecisionType>
void
FiniteStrainTensorResampleImageFilter<TImageType, TInterpolatorPrecisionType>
::RotateInterpolatedModel(const InputPixelType &interpolatedModel, vnl_matrix <double> &modelRotationMatrix,
                          InputPixelType &rotatedModel, itk::ThreadIdType threadId)
{
    anima::GetTensorFromVectorRepresentation(interpolatedModel,m_WorkMats[threadId],m_TensorDimension,true);

    anima::RotateSymmetricMatrix(m_WorkMats[threadId],modelRotationMatrix,m_TmpTensors[threadId]);

    anima::GetVectorRepresentation(m_TmpTensors[threadId],rotatedModel,m_VectorSize,true);
}

} // end of namespace anima
