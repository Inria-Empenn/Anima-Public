#pragma once
#include "animaFuzzyCMeansFilter.h"

#include <iostream>
#include <math.h>
#include <float.h>

#include <animaLogExpMapsUnitSphere.h>

namespace anima
{

template <class ScalarType>
FuzzyCMeansFilter <ScalarType>
::FuzzyCMeansFilter()
{
    m_ClassesMembership.clear();
    m_Centroids.clear();
    m_InputData.clear();

    m_NbClass = 0;
    m_NbInputs = 0;
    m_NDim = 0;
    m_MaxIterations = 100;

    m_Verbose = true;
    m_SpectralClusterInit = false;
    m_SphericalAverageType = Euclidean;

    m_RelStopCriterion = 1.0e-4;
    m_MValue = 2;
}

template <class ScalarType>
FuzzyCMeansFilter <ScalarType>
::~FuzzyCMeansFilter()
{
    m_ClassesMembership.clear();
    m_Centroids.clear();
    m_InputData.clear();
}

template <class ScalarType>
void
FuzzyCMeansFilter <ScalarType>
::SetInputData(DataHolderType &data)
{
    if (data.size() == 0)
        return;

    m_InputData = data;
    m_NbInputs = m_InputData.size();
    m_NDim = m_InputData[0].size();
}

template <class ScalarType>
void
FuzzyCMeansFilter <ScalarType>
::Update()
{
    if (m_NbClass > m_NbInputs)
        throw itk::ExceptionObject(__FILE__,__LINE__,"More classes than inputs...",ITK_LOCATION);

    if (m_DataWeights.size() != m_NbInputs)
    {
        m_DataWeights.resize(m_NbInputs);
        std::fill(m_DataWeights.begin(),m_DataWeights.end(),1.0 / m_NbInputs);
    }

    InitializeCMeansFromData();
    DataHolderType oldMemberships = m_ClassesMembership;

    unsigned int itncount = 0;
    bool continueLoop = true;

    while ((itncount < m_MaxIterations)&&(continueLoop))
    {
        itncount++;

        if (m_Verbose)
            std::cout << "Iteration " << itncount << "..." << std::endl;

        ComputeCentroids();
        UpdateMemberships();

        continueLoop = !endConditionReached(oldMemberships);
        oldMemberships = m_ClassesMembership;
    }
}

template <class ScalarType>
void
FuzzyCMeansFilter <ScalarType>
::ComputeCentroids()
{
    if (m_PowMemberships.size() != m_NbInputs)
        m_PowMemberships.resize(m_NbInputs);

    for (unsigned int i = 0;i < m_NbInputs;++i)
    {
        if (m_PowMemberships[i].size() != m_NbClass)
            m_PowMemberships[i].resize(m_NbClass);

        for (unsigned int j = 0;j < m_NbClass;++j)
            m_PowMemberships[i][j] = std::pow(m_ClassesMembership[i][j],m_MValue);
    }

    if (m_TmpVector.size() != m_NDim)
        m_TmpVector.resize(m_NDim);

    if (m_TmpWeights.size() != m_NbInputs)
        m_TmpWeights.resize(m_NbInputs);

    for (unsigned int i = 0;i < m_NbClass;++i)
    {
        double sumPowMemberShips = 0;
        std::fill(m_TmpVector.begin(),m_TmpVector.end(),0.0);

        for (unsigned int j = 0;j < m_NbInputs;++j)
        {
            sumPowMemberShips += m_DataWeights[j] * m_PowMemberships[j][i];
            for (unsigned int k = 0;k < m_NDim;++k)
                m_TmpVector[k] += m_DataWeights[j] * m_PowMemberships[j][i] * m_InputData[j][k];
        }

        for (unsigned int k = 0;k < m_NDim;++k)
            m_TmpVector[k] /= sumPowMemberShips;

        switch (m_SphericalAverageType)
        {
            case Euclidean:
                m_Centroids[i] = m_TmpVector;
                break;

            case ApproximateSpherical:
            {
                double tmpSum = 0;
                for (unsigned int k = 0;k < m_NDim;++k)
                    tmpSum += m_TmpVector[k] * m_TmpVector[k];

                tmpSum = std::sqrt(tmpSum);
                for (unsigned int k = 0;k < m_NDim;++k)
                    m_TmpVector[k] /= tmpSum;

                m_Centroids[i] = m_TmpVector;
                break;
            }

            case Spherical:
            default:
            {
                double tmpSum = 0;
                for (unsigned int k = 0;k < m_NDim;++k)
                    tmpSum += m_TmpVector[k] * m_TmpVector[k];

                tmpSum = std::sqrt(tmpSum);
                for (unsigned int k = 0;k < m_NDim;++k)
                    m_TmpVector[k] /= tmpSum;

                for (unsigned int k = 0;k < m_NbInputs;++k)
                    m_TmpWeights[k] = m_DataWeights[k] * m_PowMemberships[k][i];

                anima::ComputeSphericalCentroid(m_InputData,m_Centroids[i],m_TmpVector,m_TmpWeights,&m_WorkLogVector,&m_WorkVector);

                break;
            }
        }
    }
}

template <class ScalarType>
void
FuzzyCMeansFilter <ScalarType>
::UpdateMemberships()
{
    long double powFactor = 1.0/(m_MValue - 1.0);

    if (m_DistancesPointsCentroids.size() != m_NbClass)
        m_DistancesPointsCentroids.resize(m_NbClass);

    for (unsigned int i = 0;i < m_NbInputs;++i)
    {
        unsigned int minClassIndex = 0;
        bool nullDistance = false;
        for (unsigned int j = 0;j < m_NbClass;++j)
        {
            m_DistancesPointsCentroids[j] = computeDistance(m_InputData[i],m_Centroids[j]);

            if (m_DistancesPointsCentroids[j] <= 0)
            {
                nullDistance = true;
                minClassIndex = j;
                break;
            }
        }

        if (nullDistance)
        {
            for (unsigned int j = 0;j < m_NbClass;++j)
                m_ClassesMembership[i][j] = 0;

            m_ClassesMembership[i][minClassIndex] = 1.0;
            continue;
        }

        for (unsigned int j = 0;j < m_NbClass;++j)
        {
            long double tmpVal = 0;

            if (m_MValue == 2.0)
            {
                for (unsigned int k = 0;k < m_NbClass;++k)
                    tmpVal += m_DistancesPointsCentroids[j] / m_DistancesPointsCentroids[k];
            }
            else
            {
                for (unsigned int k = 0;k < m_NbClass;++k)
                    tmpVal += std::pow(m_DistancesPointsCentroids[j] / m_DistancesPointsCentroids[k],powFactor);
            }

            m_ClassesMembership[i][j] = 1.0 / tmpVal;
        }
    }
}

template <class ScalarType>
bool
FuzzyCMeansFilter <ScalarType>
::endConditionReached(DataHolderType &oldMemberships)
{
    double absDiff = 0;

    for (unsigned int i = 0;i < m_NbInputs;++i)
    {
        for (unsigned int j = 0;j < m_NbClass;++j)
        {
            double testValue = std::abs(oldMemberships[i][j] - m_ClassesMembership[i][j]);
            if (testValue > absDiff)
                absDiff = testValue;
        }
    }

    if (absDiff > m_RelStopCriterion)
        return false;
    else
        return true;
}

template <class ScalarType>
void
FuzzyCMeansFilter <ScalarType>
::InitializeCMeansFromData()
{
    m_Centroids.resize(m_NbClass);

    if (!m_SpectralClusterInit)
    {
        for (unsigned int i = 0;i < m_NbClass;++i)
            m_Centroids[i] = m_InputData[i];

        if (m_ClassesMembership.size() != m_NbInputs)
        {
            m_ClassesMembership.resize(m_NbInputs);

            double fixVal = 0.95;
            VectorType tmpVec(m_NbClass,0);
            for (unsigned int i = 0;i < m_NbInputs;++i)
            {
                unsigned int tmp = i % m_NbClass;
                m_ClassesMembership[i] = tmpVec;
                m_ClassesMembership[i][tmp] = fixVal;
                for (unsigned j = 0;j < m_NbClass;++j)
                {
                    if (j != tmp)
                        m_ClassesMembership[i][j] = (1.0 - fixVal)/(m_NbClass - 1.0);
                }
            }
        }
    }
    else
    {
        m_Centroids[0] = m_InputData[0];
        std::vector <unsigned int> alreadyIn(m_NbClass,0);

        for (unsigned int i = 1;i < m_NbClass;++i)
        {
            double minCrossProd = DBL_MAX;
            unsigned int minIndex = 0;
            for (unsigned int j = 0;j < m_NbInputs;++j)
            {
                bool useIt = true;
                for (unsigned int k = 0;k < i;++k)
                {
                    if (alreadyIn[k] == j)
                    {
                        useIt = false;
                        break;
                    }
                }

                if (useIt)
                {
                    double maxCrossProd = 0;
                    for (unsigned int l = 0;l < i;++l)
                    {
                        double crossProd = 0;
                        for (unsigned int k = 0;k < m_NDim;++k)
                            crossProd += m_InputData[j][k]*m_Centroids[l][k];

                        if (crossProd > maxCrossProd)
                            maxCrossProd = crossProd;
                    }

                    if (maxCrossProd < minCrossProd)
                    {
                        minCrossProd = maxCrossProd;
                        minIndex = j;
                    }
                }
            }

            m_Centroids[i] = m_InputData[minIndex];
            alreadyIn[i] = minIndex;
        }

        //Centroids initialized, now compute memberships
        if (m_ClassesMembership.size() != m_NbInputs)
        {
            m_ClassesMembership.resize(m_NbInputs);
            VectorType tmpVec(m_NbClass,0);
            for (unsigned int i = 0;i < m_NbInputs;++i)
                m_ClassesMembership[i] = tmpVec;
        }
        this->UpdateMemberships();
    }
}

template <class ScalarType>
void
FuzzyCMeansFilter <ScalarType>
::InitializeClassesMemberships(DataHolderType &classM)
{
    if (classM.size() == m_NbInputs)
    {
        m_ClassesMembership.resize(m_NbInputs);

        for (unsigned int i = 0;i < m_NbInputs;++i)
            m_ClassesMembership[i] = classM[i];
    }
}

template <class ScalarType>
long double
FuzzyCMeansFilter <ScalarType>
::computeDistance(VectorType &vec1, VectorType &vec2)
{
    long double resVal = 0;

    if (m_SphericalAverageType != Euclidean)
    {
        long double dotProd = 0;
        for (unsigned int i = 0;i < m_NDim;++i)
            dotProd += vec1[i]*vec2[i];

        if (dotProd > 1)
            dotProd = 1;

        resVal = std::abs(std::acos(dotProd));
    }
    else
    {
        resVal = 0;
        for (unsigned int i = 0;i < m_NDim;++i)
            resVal += (vec1[i] - vec2[i])*(vec1[i] - vec2[i]);
    }

    return resVal;
}

} // end namespace anima
