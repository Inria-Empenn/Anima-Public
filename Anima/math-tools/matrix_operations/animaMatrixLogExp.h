#pragma once

#include <vector>
#include <vnl/vnl_matrix.h>
#include <math.h>

#include <itkTransform.h>
#include <itkMatrixOffsetTransformBase.h>
#include <itkMultiThreaderBase.h>

namespace anima
{
/* Implementation of square matrix logarithm and exponential. First implemented by V. Arsigny. */

//! Gets the square root of matrix m
template <class T> vnl_matrix <T> GetSquareRoot(const vnl_matrix <T> & m, const double precision, vnl_matrix <T> & resultM);

//! Final part of the computation of the log. Estimates the log with a Pade approximation for a matrix m such that \|m-Id\| <= 0.5.
template <class T> vnl_matrix <T> GetPadeLogarithm(const vnl_matrix <double> & m,const int numApprox);

//! Computation of the matrix logarithm. Algo: inverse scaling and squaring, variant proposed by Cheng et al., SIAM Matrix Anal., 2001.
template <class T> vnl_matrix <T> GetLogarithm(const vnl_matrix <T> & m, const double precision=0.00000000001, const int numApprox=1);

//! Computation of the matrix exponential. Algo: classical scaling and squaring, as in Matlab. See Higham, SIAM Matr. Anal., 2004.
template <class T> vnl_matrix <T> GetExponential(const vnl_matrix <T> & m, const int numApprox=1);

//! Computation of a 3D rotation matrix logarithm. Rodrigues' explicit formula
template <class T> void Get3DRotationLogarithm(const vnl_matrix <T> &rotationMatrix, std::vector <T> &outputAngles);

//! Computation of a 3D rotation matrix exponential. Rodrigues' explicit formula
template <class T> void Get3DRotationExponential(const std::vector <T> &angles, vnl_matrix <T> &outputRotation);

//! Class to compute many log-vectors in a multi-threaded way
template <class TInputScalarType, class TOutputScalarType, unsigned int NDimensions, unsigned int NDegreesOfFreedom>
class MatrixLoggerFilter
{
public:
    typedef MatrixLoggerFilter Self;

    typedef itk::Vector <TOutputScalarType, NDegreesOfFreedom> OutputVectorType;
    typedef typename std::vector <OutputVectorType> OutputType;

    typedef itk::Transform<TInputScalarType,NDimensions,NDimensions> BaseInputTransformType;
    typedef typename std::vector <typename BaseInputTransformType::Pointer> InputType;
    typedef itk::MatrixOffsetTransformBase<TInputScalarType,NDimensions,NDimensions> BaseInputMatrixTransformType;

    MatrixLoggerFilter()
    {
        m_Modified = true;
        m_InputTransforms.clear();

        m_UseRigidTransforms = false;

        m_NumberOfThreads = itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads();
    }

    ~MatrixLoggerFilter() {}

    void SetInput(InputType &input) {m_InputTransforms = input;m_Modified = true;}
    void SetNumberOfWorkUnits(unsigned int &numThreads) {m_NumberOfThreads = numThreads;}
    void SetUseRigidTransforms(bool rigTrsfs)
    {
        if (m_UseRigidTransforms != rigTrsfs)
            m_Modified = true;

        m_UseRigidTransforms = rigTrsfs;
    }

    void Update();

    OutputType &GetOutput()
    {
        if (m_Modified)
            this->Update();

        return m_Output;
    }

protected:
    struct ThreadedLogData
    {
        Self* loggerFilter;
    };

    static ITK_THREAD_RETURN_TYPE ThreadedLogging(void *arg);
    void InternalLogger(unsigned int threadId, unsigned int numThreads);

private:
    InputType m_InputTransforms;
    OutputType m_Output;

    unsigned int m_NumberOfThreads;
    bool m_UseRigidTransforms;
    bool m_Modified;
};

} // end of namespace anima

#include "animaMatrixLogExp.hxx"
