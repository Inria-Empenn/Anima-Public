#pragma once

#include <vector>

namespace anima
{

/** \class FuzzyCMeansFilter
 * Provides an implementation of fuzzy c-means, as proposed in
 * J. C. Bezdek (1981): "Pattern Recognition with Fuzzy Objective Function Algoritms", Plenum Press, New York.
 * It contains flags for interfacing it with spectral clustering, and to compute spherical distances rather than
 * Euclidean distances.
 */
template <class ScalarType>
class FuzzyCMeansFilter
{
public:
    typedef std::vector <ScalarType> VectorType;
    typedef std::vector <VectorType> DataHolderType;

    enum CentroidAverageType
    {
        Euclidean = 0,
        ApproximateSpherical,
        Spherical
    };

    FuzzyCMeansFilter();
    virtual ~FuzzyCMeansFilter();

    void SetInputData(DataHolderType &data);
    void SetDataWeights(VectorType &val) {m_DataWeights = val;}
    void SetNbClass(unsigned int nbC) {m_NbClass = nbC;}
    void SetFlagSpectralClustering(bool flag) {m_SpectralClusterInit = flag;}
    void SetMaxIterations(unsigned int mIt) {m_MaxIterations = mIt;}
    void SetRelStopCriterion(double rC) {m_RelStopCriterion = rC;}
    void SetMValue(double mV) {m_MValue = mV;}
    void SetSphericalAverageType(CentroidAverageType spher) {m_SphericalAverageType = spher;}

    void ComputeCentroids();
    void UpdateMemberships();

    void InitializeCMeansFromData();
    void InitializeClassesMemberships(DataHolderType &classM);
    void ResetClassesMemberships() {m_ClassesMembership.clear();}

    bool endConditionReached(DataHolderType &oldMemberships);
    void SetVerbose(bool verb) {m_Verbose = verb;}

    void Update();

    VectorType &GetCentroid(unsigned int i) {return m_Centroids[i];}
    VectorType &GetClassesMembership(unsigned int i) {return m_ClassesMembership[i];}

private:
    long double computeDistance(VectorType &vec1, VectorType &vec2);

    DataHolderType m_ClassesMembership;
    DataHolderType m_Centroids;
    DataHolderType m_InputData;
    VectorType m_DataWeights;

    unsigned int m_NbClass, m_NbInputs;
    unsigned int m_NDim;
    unsigned int m_MaxIterations;

    bool m_Verbose;
    bool m_SpectralClusterInit;

    CentroidAverageType m_SphericalAverageType;

    double m_RelStopCriterion;
    double m_MValue;

    // Internal work values
    std::vector <long double> m_DistancesPointsCentroids;
    DataHolderType m_PowMemberships;
    VectorType m_TmpVector;
    VectorType m_TmpWeights;

    // Internal variable for faster pseudo-spherical average computation
    VectorType m_WorkVector, m_WorkLogVector;
};

} // end namespace anima

#include "animaFuzzyCMeansFilter.hxx"
