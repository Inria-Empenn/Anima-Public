#pragma once

#include <animaBaseTensorTools.h>

namespace anima
{

template <typename TInputScalarType, unsigned int Dimension, typename TInterpolatorPrecisionType>
void
FiniteStrainTensorResampleImageFilter<TInputScalarType, Dimension, TInterpolatorPrecisionType>
::RotateInterpolatedModel(const InputPixelType &interpolatedModel, vnl_matrix <double> &modelRotationMatrix,
                          InputPixelType &rotatedModel)
{
    vnl_matrix <double> workMat(m_TensorDimension,m_TensorDimension);
    vnl_matrix <double> tmpTensor(m_TensorDimension,m_TensorDimension,0);

    anima::GetTensorFromVectorRepresentation(interpolatedModel,workMat,m_TensorDimension,true);

    anima::RotateSymmetricMatrix(workMat,modelRotationMatrix,tmpTensor);

    anima::GetVectorRepresentation(tmpTensor,rotatedModel,m_VectorSize,true);
}

} // end of namespace anima
