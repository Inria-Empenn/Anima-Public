#pragma once
#include "animaKMeansFilter.h"

namespace anima {

template <class DataType, unsigned int PointDimension>
KMeansFilter <DataType,PointDimension>::
KMeansFilter()
{
    m_ClassesMembership.clear();
    m_Centroids.clear();
    m_InputData.clear();
    m_NumberPerClass.clear();

    m_NbClass = 0;
    m_NbInputs = 0;
    m_MaxIterations = 100;

    m_Verbose = true;
}

template <class DataType, unsigned int PointDimension>
KMeansFilter <DataType,PointDimension>::
~KMeansFilter()
{
}

template <class DataType, unsigned int PointDimension>
void
KMeansFilter <DataType,PointDimension>::
SetInputData(DataHolderType &data)
{
    if (data.size() == 0)
        return;

    m_InputData = data;
    m_NbInputs = m_InputData.size();
}

template <class DataType, unsigned int PointDimension>
void
KMeansFilter <DataType,PointDimension>::
Update()
{
    if (m_NbClass > m_NbInputs)
        throw itk::ExceptionObject(__FILE__, __LINE__,"More classes than inputs...",ITK_LOCATION);

    this->InitializeKMeansFromData();
    MembershipType oldMemberships = m_ClassesMembership;
    unsigned int itncount = 0;
    bool continueLoop = true;

    while ((itncount < m_MaxIterations)&&(continueLoop))
    {
        itncount++;

        if (m_Verbose)
            std::cout << "Iteration " << itncount << "..." << std::endl;

        this->ComputeCentroids();
        this->UpdateMemberships();

        continueLoop = !this->endConditionReached(oldMemberships);
        oldMemberships = m_ClassesMembership;
    }
}

template <class DataType, unsigned int PointDimension>
void
KMeansFilter <DataType,PointDimension>::
ComputeCentroids()
{
    for (unsigned int i = 0;i < m_NbClass;++i)
        m_Centroids[i].Fill(0);

    for (unsigned int j = 0;j < m_NbInputs;++j)
    {
        for (unsigned int k = 0;k < PointDimension;++k)
            m_Centroids[m_ClassesMembership[j]][k] += m_InputData[j][k];
    }

    for (unsigned int j = 0;j < m_NbClass;++j)
    {
        if (m_NumberPerClass[j] != 0)
        {
            for (unsigned int k = 0;k < PointDimension;++k)
                m_Centroids[j][k] /= m_NumberPerClass[j];
        }
    }
}

template <class DataType, unsigned int PointDimension>
void
KMeansFilter <DataType,PointDimension>::
UpdateMemberships()
{
    std::fill(m_NumberPerClass.begin(),m_NumberPerClass.end(),0);
    for (unsigned int i = 0;i < m_NbInputs;++i)
    {
        unsigned int bestClass = 0;
        double bestDistance = this->computeDistance(m_InputData[i],m_Centroids[0]);

        for (unsigned int j = 1;j < m_NbClass;++j)
        {
            double tmpDist = this->computeDistance(m_InputData[i],m_Centroids[j]);
            if (tmpDist < bestDistance)
            {
                bestDistance = tmpDist;
                bestClass = j;
            }
        }

        m_ClassesMembership[i] = bestClass;
        ++m_NumberPerClass[bestClass];
    }
}

template <class DataType, unsigned int PointDimension>
bool
KMeansFilter <DataType,PointDimension>::
endConditionReached(MembershipType &oldMemberships)
{
    for (unsigned int i = 0;i < m_NbInputs;++i)
    {
        if (oldMemberships[i] != m_ClassesMembership[i])
            return false;
    }

    return true;
}

template <class DataType, unsigned int PointDimension>
void
KMeansFilter <DataType,PointDimension>::
InitializeKMeansFromData()
{
    m_Centroids.clear();

    for (unsigned int i = 0;i < m_NbClass;++i)
        m_Centroids.push_back(m_InputData[i]);

    //Centroids initialized, now compute memberships
    if (m_ClassesMembership.size() != m_NbInputs)
    {
        m_ClassesMembership.resize(m_NbInputs);
        std::fill(m_ClassesMembership.begin(),m_ClassesMembership.end(),0);
        
        this->UpdateMemberships();
    }
}

template <class DataType, unsigned int PointDimension>
void
KMeansFilter <DataType,PointDimension>::
InitializeClassesMemberships(MembershipType &classM)
{
    if (classM.size() == m_NbInputs)
        m_ClassesMembership = classM;

    for (unsigned int i = 0;i < m_NbInputs;++i)
        m_NumberPerClass[m_ClassesMembership[i]]++;
}

template <class DataType, unsigned int PointDimension>
double
KMeansFilter <DataType,PointDimension>::
computeDistance(VectorType &vec1, VectorType &vec2)
{
    double resVal = 0;

    for (unsigned int i = 0;i < PointDimension;++i)
        resVal += (vec1[i] - vec2[i])*(vec1[i] - vec2[i]);

    return resVal;
}

} // end namespace anima
