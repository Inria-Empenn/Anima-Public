#pragma once
#include <animaTensorResampleImageFilter.h>

#include <animaBaseTensorTools.h>

namespace anima
{

template <typename TImageType, typename TInterpolatorPrecisionType>
void
TensorResampleImageFilter<TImageType, TInterpolatorPrecisionType>
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

    if (!this->GetFiniteStrainReorientation())
    {
        m_WorkEigenValues.resize(this->GetNumberOfThreads());
        m_WorkEigenVectors.resize(this->GetNumberOfThreads());
        m_WorkPPDOrientationMatrices.resize(this->GetNumberOfThreads());
    }
}

template <typename TImageType, typename TInterpolatorPrecisionType>
void
TensorResampleImageFilter<TImageType, TInterpolatorPrecisionType>
::ReorientInterpolatedModel(const InputPixelType &interpolatedModel, vnl_matrix <double> &modelOrientationMatrix,
                            InputPixelType &rotatedModel, itk::ThreadIdType threadId)
{
    anima::GetTensorFromVectorRepresentation(interpolatedModel,m_WorkMats[threadId],m_TensorDimension,true);

    if (this->GetFiniteStrainReorientation())
        anima::RotateSymmetricMatrix(m_WorkMats[threadId],modelOrientationMatrix,m_TmpTensors[threadId]);
    else
    {
        typedef itk::Matrix <double,3,3> MatrixType;
        typedef vnl_vector_fixed <double,3> VectorType;
        itk::SymmetricEigenAnalysis < MatrixType, VectorType, MatrixType> eigenComputer(3);
        eigenComputer.ComputeEigenValuesAndVectors(m_WorkMats[threadId],m_WorkEigenValues[threadId],
                                                   m_WorkEigenVectors[threadId]);

        anima::ExtractPPDRotationFromJacobianMatrix(modelOrientationMatrix,m_WorkPPDOrientationMatrices[threadId],m_WorkEigenVectors[threadId]);
        anima::RotateSymmetricMatrix(m_WorkMats[threadId],m_WorkPPDOrientationMatrices[threadId],m_TmpTensors[threadId]);
    }

    anima::GetVectorRepresentation(m_TmpTensors[threadId],rotatedModel,m_VectorSize,true);
}

} // end of namespace anima
