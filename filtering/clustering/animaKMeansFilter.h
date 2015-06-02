#pragma once

#include <vector>

namespace anima {

template <class DataType, unsigned int PointDimension>
class KMeansFilter
{
public:
    typedef DataType VectorType;
    typedef std::vector < VectorType > DataHolderType;
    typedef std::vector < unsigned int > MembershipType;

    KMeansFilter();
    virtual ~KMeansFilter();

    void SetInputData(DataHolderType &data);
    void SetNumberOfClasses(unsigned int nbC)
    {
        m_NbClass = nbC;
        m_NumberPerClass.resize(m_NbClass);
    }

    void SetMaxIterations(unsigned int mIt) {m_MaxIterations = mIt;}

    void ComputeCentroids();
    void UpdateMemberships();

    void InitializeKMeansFromData();
    void InitializeClassesMemberships(MembershipType &classM);
    void ResetClassesMemberships() {m_ClassesMembership.clear();}

    bool endConditionReached(MembershipType &oldMemberships);
    void SetVerbose(bool verb) {m_Verbose = verb;}

    void Update();

    VectorType GetCentroid(unsigned int i) {return m_Centroids[i];}
    unsigned int GetClassMembership(unsigned int i) {return m_ClassesMembership[i];}
    MembershipType &GetClassesMemberships() {return m_ClassesMembership;}

    unsigned int GetNumberPerClass(unsigned int i) {return m_NumberPerClass[i];}

private:
    double computeDistance(VectorType &vec1, VectorType &vec2);

    MembershipType m_ClassesMembership;
    DataHolderType m_Centroids;
    DataHolderType m_InputData;

    std::vector <unsigned int> m_NumberPerClass;

    unsigned int m_NbClass, m_NbInputs;
    unsigned int m_MaxIterations;

    bool m_Verbose;
};

} // end namespace anima

#include "animaKMeansFilter.hxx"
