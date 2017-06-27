#pragma once
#include "animaSpectralClusteringFilter.h"

#include <itkSymmetricEigenAnalysis.h>
#include <vnl/vnl_diag_matrix.h>

namespace anima
{

template <class ScalarType>
SpectralClusteringFilter <ScalarType>
::SpectralClusteringFilter()
{
    m_ClassesMembership.clear();
    m_Centroids.clear();
    m_SpectralVectors.clear();
    m_InputData.set_size(1,1);

    m_NbClass = 0;
    m_MaxIterations = 100;

    m_CMeansAverageType = CMeansFilterType::ApproximateSpherical;

    m_Verbose = true;
    m_RelStopCriterion = 1.0e-4;
    m_MValue = 2;

    m_SigmaWeighting = 1;
}

template <class ScalarType>
void
SpectralClusteringFilter <ScalarType>
::ResetOutputs()
{
    m_ClassesMembership.clear();
    m_Centroids.clear();
    m_SpectralVectors.clear();
}

template <class ScalarType>
SpectralClusteringFilter <ScalarType>
::~SpectralClusteringFilter()
{
    this->ResetOutputs();
}

template <class ScalarType>
void
SpectralClusteringFilter <ScalarType>
::SetInputData(MatrixType &data)
{
    if (data.rows() == 0)
        return;

    m_InputData = data;
}

template <class ScalarType>
void
SpectralClusteringFilter <ScalarType>
::Update()
{
    if (m_NbClass > m_InputData.rows())
        throw itk::ExceptionObject(__FILE__,__LINE__,"More classes than inputs...",ITK_LOCATION);

    if (m_DataWeights.size() != m_InputData.rows())
    {
        m_DataWeights.resize(m_InputData.rows());
        std::fill(m_DataWeights.begin(),m_DataWeights.end(),1.0 / m_InputData.rows());
    }

    this->ComputeSpectralVectors();

    CMeansFilterType mainFilter;
    mainFilter.SetNbClass(m_NbClass);
    mainFilter.SetMaxIterations(m_MaxIterations);
    mainFilter.SetRelStopCriterion(m_RelStopCriterion);
    mainFilter.SetMValue(m_MValue);

    mainFilter.SetInputData(m_SpectralVectors);
    mainFilter.SetDataWeights(m_DataWeights);
    mainFilter.SetVerbose(m_Verbose);
    mainFilter.SetFlagSpectralClustering(true);
    mainFilter.SetSphericalAverageType(m_CMeansAverageType);

    mainFilter.Update();

    unsigned int inputSize = m_InputData.rows();
    m_ClassesMembership.clear();

    for (unsigned int i = 0;i < inputSize;++i)
        m_ClassesMembership.push_back(mainFilter.GetClassesMembership(i));

    m_Centroids.clear();
    for (unsigned int i = 0;i < m_NbClass;++i)
        m_Centroids.push_back(mainFilter.GetCentroid(i));
}

template <class ScalarType>
void
SpectralClusteringFilter <ScalarType>
::InitializeSigmaFromDistances()
{
    double sigmaTmp = 0;

    unsigned int inputSize = m_InputData.rows();

    unsigned int numPts = inputSize*(inputSize + 1)/2 - inputSize - 1;

    for (unsigned int i = 0;i < inputSize;++i)
        for (unsigned int j = i+1;j < inputSize;++j)
            sigmaTmp += m_InputData(i,j);

    if (sigmaTmp > 0)
        m_SigmaWeighting = std::sqrt(sigmaTmp/numPts);
    else
        m_SigmaWeighting = 1;
}

template <class ScalarType>
void
SpectralClusteringFilter <ScalarType>
::ComputeSpectralVectors()
{
    unsigned int inputSize = m_InputData.rows();
    MatrixType W(inputSize,inputSize,0);
    std::vector <double> dValues(inputSize,0);

    for (unsigned int i = 0;i < inputSize;++i)
        for (unsigned int j = i+1;j < inputSize;++j)
        {
            W(i,j) = std::exp(- m_InputData(i,j) / (2.0 * m_SigmaWeighting * m_SigmaWeighting));
            W(j,i) = W(i,j);
        }

    for (unsigned int i = 0;i < inputSize;++i)
    {
        dValues[i] = 0;
        for (unsigned int j = 0;j < inputSize;++j)
            dValues[i] += W(i,j);

        dValues[i] = 1.0/std::sqrt(dValues[i]);
    }

    for (unsigned int i = 0;i < inputSize;++i)
        for (unsigned int j = i+1;j < inputSize;++j)
        {
            W(i,j) *= dValues[i]*dValues[j];
            W(j,i) = W(i,j);
        }

    // Compute eigenvalues and vectors, keep only m_NbClasses largest values and corresponding truncated vectors
    vnl_diag_matrix<ScalarType> eigVals(inputSize);
    eigVals.fill(0.0);
    MatrixType eigVecs(inputSize,inputSize);
    eigVecs.fill(0.0);

    typedef itk::SymmetricEigenAnalysis <MatrixType, vnl_diag_matrix<ScalarType>, MatrixType> EigenAnalysisType;
    EigenAnalysisType eigen(inputSize);
    eigen.ComputeEigenValuesAndVectors(W,eigVals,eigVecs);

    m_SpectralVectors.resize(inputSize);
    VectorType tmpVec(m_NbClass);
    for (unsigned int i = 0;i < inputSize;++i)
    {
        for (unsigned int j = 0;j < m_NbClass;++j)
            tmpVec[j] = eigVecs(inputSize - j - 1,i);

        double tmpSum = 0;
        for (unsigned int j = 0;j < m_NbClass;++j)
            tmpSum += tmpVec[j]*tmpVec[j];

        tmpSum = std::sqrt(tmpSum);
        for (unsigned int j = 0;j < m_NbClass;++j)
            tmpVec[j] /= tmpSum;

        m_SpectralVectors[i] = tmpVec;
    }
}

template <class ScalarType>
std::vector <unsigned int>
SpectralClusteringFilter <ScalarType>
::GetClassMembers(unsigned int i)
{
    std::vector <unsigned int> resVal;

    if (m_ClassesMembership.size() != m_InputData.rows())
        return resVal;

    unsigned int inputSize = m_InputData.rows();

    for (unsigned int j = 0;j < inputSize;++j)
    {
        bool isClassIMax = true;
        for (unsigned int k = 0;k < m_NbClass;++k)
        {
            if (m_ClassesMembership[j][k] > m_ClassesMembership[j][i])
            {
                isClassIMax = false;
                break;
            }
        }

        if (isClassIMax)
            resVal.push_back(j);
    }

    return resVal;
}

template <class ScalarType>
double
SpectralClusteringFilter <ScalarType>
::ComputeClustersSpreading()
{
    if (m_ClassesMembership.size() != m_InputData.rows())
        return -1;

    unsigned int inputSize = m_InputData.rows();
    std::vector <double> classDistance(m_NbClass,0);

    for (unsigned int i = 0;i < m_NbClass;++i)
    {
        double sumWeights = 0;
        for (unsigned int j = 0;j < inputSize;++j)
        {
            sumWeights += m_ClassesMembership[j][i];

            for (unsigned int k = 0;k < m_NbClass;++k)
                classDistance[i] += m_ClassesMembership[j][i] * std::abs(m_Centroids[i][k] - m_SpectralVectors[j][k]);
        }

        classDistance[i] /= sumWeights;
    }

    double resVal = 0;
    for (unsigned int i = 0;i < m_NbClass;++i)
        resVal += classDistance[i];

    return resVal/m_NbClass;
}

} // end namespace anima
