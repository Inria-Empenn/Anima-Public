#pragma once

#include "animaModularityClusteringFilter.h"

#include <itkSymmetricEigenAnalysis.h>

namespace anima {

template <class DataType>
const double ModularityClusteringFilter <DataType>::m_ZeroThreshold = 1.0e-6;

template <class DataType>
ModularityClusteringFilter <DataType>::
ModularityClusteringFilter()
{
    m_ClassMemberships.clear();
    m_ReverseClassMemberships.clear();
    m_ConnectivityMatrix.clear();
    m_ModularityMatrix.clear();
    m_SubModularityMatrix.clear();
    m_DegreeMatrix.clear();
    m_SubDegreeMatrix.clear();
    m_PreventClusterFromSplitting.clear();
    m_SubMembership1.clear();
    m_SubMembership2.clear();
    m_Membership.clear();

    m_NbInputs = 0;
    m_NumberOfClusters = 0;

    m_TotalNumberOfEdges = 0.0;
}

template <class DataType>
ModularityClusteringFilter <DataType>::
~ModularityClusteringFilter()
{
}

template <class DataType>
void
ModularityClusteringFilter <DataType>::
SetInputData(DataHolderType &data)
{
    if (data.rows() == 0)
        return;

    m_ConnectivityMatrix = data;

    m_NbInputs = m_ConnectivityMatrix.rows();

    m_TotalNumberOfEdges = 0.0;
    m_DegreeMatrix.set_size(m_NbInputs);

    for (unsigned int i = 0;i < m_NbInputs;++i)
    {
        double t = m_ConnectivityMatrix.get_row(i).sum();
        m_DegreeMatrix(i) = t;
        m_TotalNumberOfEdges += t;
    }

    m_TotalNumberOfEdges /= 2.0;
}

template <class DataType>
void
ModularityClusteringFilter <DataType>::
InitializeClassMemberships()
{
    m_NumberOfClusters = 1;

    m_ClassMemberships.resize(m_NbInputs);
    m_ReverseClassMemberships.resize(m_NumberOfClusters);
    m_ReverseClassMemberships[0].resize(m_NbInputs);
    for (unsigned int i = 0;i < m_NbInputs;++i)
    {
        m_ClassMemberships[i] = 0;
        m_ReverseClassMemberships[0][i] = i;
    }

    m_PreventClusterFromSplitting.resize(m_NumberOfClusters);
    m_PreventClusterFromSplitting[0] = false;
}

template <class DataType>
void
ModularityClusteringFilter <DataType>::
GetModularityMatrix()
{
    m_ModularityMatrix.set_size(m_NbInputs,m_NbInputs);

    for (unsigned int i = 0;i < m_NbInputs;++i)
    {
        for (unsigned int j = i;j < m_NbInputs;++j)
        {
            double t = m_ConnectivityMatrix(i,j) - m_DegreeMatrix(i) * m_DegreeMatrix(j) / (2.0 * m_TotalNumberOfEdges);
            m_ModularityMatrix(i,j) = t;

            if (i != j)
                m_ModularityMatrix(j,i) = t;
        }
    }
}

template <class DataType>
void
ModularityClusteringFilter <DataType>::
GetSubModularityMatrix(unsigned int parentIndex)
{
    m_Membership = m_ReverseClassMemberships[parentIndex];
    unsigned int numData = m_Membership.size();

    m_SubDegreeMatrix.set_size(numData);
    m_SubDegreeMatrix.fill(0.0);

    for (unsigned int i = 0;i < numData;++i)
        for (unsigned int j = 0;j < numData;++j)
            m_SubDegreeMatrix(i) += m_ModularityMatrix(m_Membership[i],m_Membership[j]);

    m_SubModularityMatrix.set_size(numData,numData);

    for (unsigned int i = 0;i < numData;++i)
    {
        for (unsigned int j = i;j < numData;++j)
        {
            double t = m_ModularityMatrix(m_Membership[i],m_Membership[j]);

            m_SubModularityMatrix(i,j) = t;

            if (i == j)
                m_SubModularityMatrix(i,j) -= m_SubDegreeMatrix(i);
            else
                m_SubModularityMatrix(j,i) = t;
        }
    }
}

template <class DataType>
double
ModularityClusteringFilter <DataType>::
SpectralClustering(unsigned int parentIndex)
{
    m_Membership = m_ReverseClassMemberships[parentIndex];
    unsigned int numData = m_Membership.size();

    // Data for eigen analysis
    typedef itk::SymmetricEigenAnalysis <MatrixType,VectorType,MatrixType> eigenAnalysis;
    MatrixType eigVecs(numData,numData);
    VectorType eigVals(numData);

    this->GetSubModularityMatrix(parentIndex);

    eigenAnalysis eigSystem(numData);
    eigSystem.SetOrderEigenValues(true);
    eigSystem.ComputeEigenValuesAndVectors(m_SubModularityMatrix, eigVals, eigVecs);

    VectorType tmpVec = eigVecs.get_row(numData-1);

    for (unsigned int i = 0;i < numData;++i)
    {
        unsigned int currentParticle = m_Membership[i];

        if (tmpVec(i) < 0)
            m_SubMembership1.push_back(currentParticle);
        else
            m_SubMembership2.push_back(currentParticle);
    }

    // Compute modularity for splitting this cluster into two subclusters
    double s_i, s_j, Q = 0;
    for (unsigned int i = 0;i < numData;++i)
    {
        if (tmpVec(i) < 0)
            s_i = -1;
        else
            s_i = 1;

        for (unsigned int j = 0;j < numData;++j)
        {
            if (tmpVec(j) < 0)
                s_j = -1;
            else
                s_j = 1;

            Q += m_SubModularityMatrix(i,j) * s_i * s_j;
        }
    }

    Q /= (4.0 * m_TotalNumberOfEdges);
    return Q;
}

template <class DataType>
void
ModularityClusteringFilter <DataType>::
SplitCluster(unsigned int parentIndex)
{
    m_Membership = m_ReverseClassMemberships[parentIndex];
    unsigned int numData = m_Membership.size();

    double Q = 0;

    // perform traditional spectral clustering
    m_SubMembership1.clear();
    m_SubMembership2.clear();
    Q = this->SpectralClustering(parentIndex);

    if (Q > m_ZeroThreshold && m_SubMembership1.size() != numData && m_SubMembership2.size() != numData)
    {
        // Update class memberships
        for (unsigned int i = 0;i < m_SubMembership2.size();++i)
            m_ClassMemberships[m_SubMembership2[i]] = m_NumberOfClusters;

        m_ReverseClassMemberships[parentIndex] = m_SubMembership1;
        m_ReverseClassMemberships.push_back(m_SubMembership2);
        ++m_NumberOfClusters;

        m_PreventClusterFromSplitting.push_back(false);
    }
    else
        m_PreventClusterFromSplitting[parentIndex] = true;
}

template <class DataType>
void
ModularityClusteringFilter <DataType>::
Update()
{
    if (m_NbInputs == 0)
    {
        std::cerr << "Input data is empty." << std::endl;
        return;
    }

    this->InitializeClassMemberships();
    this->GetModularityMatrix();

    unsigned int oldNumberOfClusters = 0;
    while (m_NumberOfClusters > oldNumberOfClusters)
    {
        oldNumberOfClusters = m_NumberOfClusters;

        for (unsigned int i = 0;i < oldNumberOfClusters;++i)
            if (!m_PreventClusterFromSplitting[i])
                this->SplitCluster(i);
    }
}
} // end namespace anima


