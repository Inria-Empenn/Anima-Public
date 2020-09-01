#pragma once
#include "animaSpectralClusteringFilter.h"

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

    m_MainFilter.SetNbClass(m_NbClass);
    m_MainFilter.SetMaxIterations(m_MaxIterations);
    m_MainFilter.SetRelStopCriterion(m_RelStopCriterion);
    m_MainFilter.SetMValue(m_MValue);

    m_MainFilter.SetInputData(m_SpectralVectors);
    m_MainFilter.SetDataWeights(m_DataWeights);
    m_MainFilter.SetVerbose(m_Verbose);
    m_MainFilter.SetFlagSpectralClustering(true);
    m_MainFilter.SetSphericalAverageType(m_CMeansAverageType);

    m_MainFilter.Update();

    unsigned int inputSize = m_InputData.rows();
    m_ClassesMembership.resize(inputSize);

    for (unsigned int i = 0;i < inputSize;++i)
        m_ClassesMembership[i] = m_MainFilter.GetClassesMembership(i);

    m_Centroids.resize(m_NbClass);
    for (unsigned int i = 0;i < m_NbClass;++i)
        m_Centroids[i] = m_MainFilter.GetCentroid(i);
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
    m_WMatrix.set_size(inputSize,inputSize);
    m_DValues.resize(inputSize);

    for (unsigned int i = 0;i < inputSize;++i)
    {
        for (unsigned int j = i+1;j < inputSize;++j)
        {
            m_WMatrix(i,j) = std::exp(- m_InputData(i,j) / (2.0 * m_SigmaWeighting * m_SigmaWeighting));
            m_WMatrix(j,i) = m_WMatrix(i,j);
        }
    }

    for (unsigned int i = 0;i < inputSize;++i)
    {
        m_DValues[i] = 0;
        for (unsigned int j = 0;j < inputSize;++j)
            m_DValues[i] += m_WMatrix(i,j);

        m_DValues[i] = 1.0/std::sqrt(m_DValues[i]);
    }

    for (unsigned int i = 0;i < inputSize;++i)
    {
        for (unsigned int j = i+1;j < inputSize;++j)
        {
            m_WMatrix(i,j) *= m_DValues[i] * m_DValues[j];
            m_WMatrix(j,i) = m_WMatrix(i,j);
        }
    }

    // Compute eigenvalues and vectors, keep only m_NbClasses largest values and corresponding truncated vectors
    m_EigenAnalyzer.SetDimension(inputSize);
    m_EigenAnalyzer.SetOrder(inputSize);
    m_EigVals.set_size(inputSize);
    m_EigVecs.set_size(inputSize,inputSize);
    m_EigenAnalyzer.ComputeEigenValuesAndVectors(m_WMatrix,m_EigVals,m_EigVecs);

    m_SpectralVectors.resize(inputSize);
    m_WorkVec.resize(m_NbClass);
    for (unsigned int i = 0;i < inputSize;++i)
    {
        for (unsigned int j = 0;j < m_NbClass;++j)
            m_WorkVec[j] = m_EigVecs.get(inputSize - j - 1,i);

        double tmpSum = 0;
        for (unsigned int j = 0;j < m_NbClass;++j)
            tmpSum += m_WorkVec[j]*m_WorkVec[j];

        tmpSum = std::sqrt(tmpSum);
        for (unsigned int j = 0;j < m_NbClass;++j)
            m_WorkVec[j] /= tmpSum;

        m_SpectralVectors[i] = m_WorkVec;
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
