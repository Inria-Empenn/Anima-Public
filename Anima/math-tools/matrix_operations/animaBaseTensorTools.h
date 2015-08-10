#pragma once

#include <itkMatrix.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>

#include <itkVariableLengthVector.h>
#include <itkMatrixOffsetTransformBase.h>

namespace anima
{

template <class T> void GetTensorLogarithm(const vnl_matrix <T> &tensor, vnl_matrix <T> &log_tensor);
template <class T> void GetTensorExponential(const vnl_matrix <T> &log_tensor, vnl_matrix <T> &tensor);
template <class T> void GetTensorPower(const vnl_matrix <T> &tensor, vnl_matrix <T> &outputTensor, double powerValue);

template <class T1, class T2> void GetVectorRepresentation(const vnl_matrix <T1> &tensor, itk::VariableLengthVector <T2> &vector,
                                                           unsigned int vecDim = 0, bool scale = false);

template <class T1, class T2> void GetTensorFromVectorRepresentation(const itk::VariableLengthVector <T1> &vector,
                                                                     vnl_matrix <T2> &tensor, unsigned int tensDim = 0,
                                                                     bool scale = false);

template <typename PixelRawType, unsigned int ImageDimension>
void ExtractRotationFromMatrixTransform(itk::MatrixOffsetTransformBase <PixelRawType,ImageDimension,ImageDimension> *trsf,
                                        vnl_matrix <double> &rotationMatrix, vnl_matrix <double> &tmpMat);

//! Recompose tensor from values extracted using SymmetricEigenAnalysis (vnl_symmetric_eigensystem transposes all this)
template <class T1, class T2>
void RecomposeTensor(vnl_diag_matrix <T1> &eigs, vnl_matrix <T1> &eigVecs, vnl_matrix <T2> &resMatrix);

template <class T1, class T2, class T3>
void RotateSymmetricMatrix(T1 &tensor, T2 &rotationMatrix, T3 &rotated_tensor, unsigned int tensorDim);

template <class T1, class T2>
void RotateSymmetricMatrix(vnl_matrix <T1> &tensor, vnl_matrix <T2> &rotationMatrix, vnl_matrix <T2> &rotated_tensor);

template <class T1, class T2, unsigned int NDim>
void RotateSymmetricMatrix(itk::Matrix <T1,NDim,NDim> &tensor, itk::Matrix <T2,NDim,NDim> &rotationMatrix,
                           itk::Matrix <T2,NDim,NDim> &rotated_tensor);

template <class T1> double ovlScore(vnl_diag_matrix <T1> &eigsX, vnl_matrix <T1> &eigVecsX,
                                    vnl_diag_matrix <T1> &eigsY, vnl_matrix <T1> &eigVecsY);


} // end of namespace anima

#include "animaBaseTensorTools.hxx"
