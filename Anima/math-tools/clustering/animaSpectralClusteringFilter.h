#pragma once

#include <vector>
#include <vnl/vnl_matrix.h>

#include <animaFuzzyCMeansFilter.h>

namespace anima
{

/** \class SpectralClusteringFilter
 * \brief Provides an implementation of spectral clustering, as proposed in
 * A.Y. Ng, M.I. Jordan and Y. Weiss. "On Spectral Clustering: Analysis and an Algorithm."
 * Advances in Neural Information Processing Systems 14. 2001
 */
template <class ScalarType>
class SpectralClusteringFilter
{
public:
    typedef std::vector <ScalarType> VectorType;
    typedef vnl_matrix <ScalarType> MatrixType;

    typedef anima::FuzzyCMeansFilter <ScalarType> CMeansFilterType;
    typedef typename CMeansFilterType::CentroidAverageType CMeansAverageType;

    SpectralClusteringFilter();
    virtual ~SpectralClusteringFilter();

    //! Input data: matrix of squared distances
    void SetInputData(MatrixType &data);
    void SetDataWeights(VectorType &val) {m_DataWeights = val;}
    void SetNbClass(unsigned int nbC) {m_NbClass = nbC;}
    void SetSigmaWeighting(double sigma) {m_SigmaWeighting = sigma;}

    // Parameters for fuzzy c-means
    void SetMaxIterations(unsigned int mIt) {m_MaxIterations = mIt;}
    void SetRelStopCriterion(double rC) {m_RelStopCriterion = rC;}
    void SetMValue(double mV) {m_MValue = mV;}

    void SetCMeansAverageType(CMeansAverageType val) {m_CMeansAverageType = val;}

    void ComputeSpectralVectors();
    void Update();

    void ResetOutputs();

    double ComputeClustersSpreading();
    void InitializeSigmaFromDistances();

    void SetVerbose(bool verb) {m_Verbose = verb;}

    VectorType &GetSpectralVector(unsigned int i) {return m_SpectralVectors[i];}
    VectorType &GetClassesMembership(unsigned int i) {return m_ClassesMembership[i];}
    VectorType &GetCentroid(unsigned int i) {return m_Centroids[i];}
    std::vector <unsigned int> GetClassMembers(unsigned int i);

private:
    std::vector <VectorType> m_ClassesMembership;
    std::vector <VectorType> m_Centroids;
    std::vector <VectorType> m_SpectralVectors;
    MatrixType m_InputData;
    VectorType m_DataWeights;

    unsigned int m_NbClass;

    CMeansAverageType m_CMeansAverageType;

    bool m_Verbose;

    unsigned int m_MaxIterations;
    double m_RelStopCriterion;
    double m_MValue;

    double m_SigmaWeighting;
};

} // end namespace anima

#include "animaSpectralClusteringFilter.hxx"
