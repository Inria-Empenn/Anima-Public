#pragma once

#include <itkMatrix.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>

#include <itkVariableLengthVector.h>
#include <itkMatrixOffsetTransformBase.h>
#include <itkSymmetricEigenAnalysis.h>

#include <vnl/vnl_diag_matrix.h>
#include <vnl/vnl_matrix.h>

namespace anima
{

template <class TScalarType>
class LogEuclideanTensorCalculator : public itk::LightObject
{
public:
    using Self = LogEuclideanTensorCalculator <TScalarType>;
    using Superclass = itk::LightObject;
    using Pointer = itk::SmartPointer<Self>;
    using ConstPointer = itk::SmartPointer<const Self>;

    /** Run-time type information (and related methods) */
    itkTypeMacro(LogEuclideanTensorCalculator, itk::LightObject)

    itkNewMacro(Self)

    using EigenAnalysisType = itk::SymmetricEigenAnalysis < vnl_matrix <double>, vnl_diag_matrix<double>, vnl_matrix <double> >;

    void GetTensorLogarithm(const vnl_matrix <TScalarType> &tensor, vnl_matrix <TScalarType> &log_tensor);
    void GetTensorExponential(const vnl_matrix <TScalarType> &tensor, vnl_matrix <TScalarType> &log_tensor);
    void GetTensorPower(const vnl_matrix <TScalarType> &tensor, vnl_matrix <TScalarType> &outputTensor, double powerValue);

private:
    EigenAnalysisType m_EigenAnalyzer;
    vnl_matrix <TScalarType> m_EigVecs;
    vnl_diag_matrix <TScalarType> m_EigVals;

};

template <class T1, class T2> void GetVectorRepresentation(const vnl_matrix <T1> &tensor, itk::VariableLengthVector <T2> &vector,
                                                           unsigned int vecDim = 0, bool scale = false);

template <class T1, class T2> void GetTensorFromVectorRepresentation(const itk::VariableLengthVector <T1> &vector,
                                                                     vnl_matrix <T2> &tensor, unsigned int tensDim = 0,
                                                                     bool scale = false);

//! Recompose tensor from values extracted using SymmetricEigenAnalysis (vnl_symmetric_eigensystem transposes all this)
template <class T1, class T2>
void RecomposeTensor(vnl_diag_matrix <T1> &eigs, vnl_matrix <T1> &eigVecs, vnl_matrix <T2> &resMatrix);

template <class T> void ProjectOnTensorSpace(const vnl_matrix <T> &matrix, vnl_matrix <T> &tensor);

template <typename RealType>
void ExtractRotationFromJacobianMatrix(vnl_matrix <RealType> &jacobianMatrix, vnl_matrix <RealType> &rotationMatrix,
                                       vnl_matrix <RealType> &tmpMat);

template <typename RealType, typename MatrixType>
void ExtractPPDRotationFromJacobianMatrix(vnl_matrix <RealType> &jacobianMatrix, vnl_matrix <RealType> &rotationMatrix,
                                          MatrixType &eigenVectors);

//! Rotates a symmetric matrix by performing R^T D R where R is a rotation matrix and D the symmetric matrix
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
