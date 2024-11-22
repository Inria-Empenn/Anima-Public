#pragma once

#include <animaNonLocalPatchBaseSearcher.h>

namespace anima
{

/**
 * Non local means patch searcher: computes actual patch matches weights and conformity.
 * Supports only one comparison image that is the same as the input
 */
template <class ImageType, class DataImageType>
class NonLocalMeansPatchSearcher : public anima::NonLocalPatchBaseSearcher <ImageType>
{
public:
    typedef typename DataImageType::Pointer DataImagePointer;
    typedef anima::NonLocalPatchBaseSearcher <ImageType> Superclass;
    typedef typename Superclass::ImageRegionType ImageRegionType;
    typedef typename Superclass::IndexType IndexType;

    NonLocalMeansPatchSearcher();
    virtual ~NonLocalMeansPatchSearcher() {}

    void SetBetaParameter(double arg) {m_BetaParameter = arg;}
    void SetNoiseCovariance(double arg) {m_NoiseCovariance = arg;}
    void SetMeanMinThreshold(double arg) {m_MeanMinThreshold = arg;}
    void SetVarMinThreshold(double arg) {m_VarMinThreshold = arg;}

    void SetMeanImage(DataImageType *arg) {m_MeanImage = arg;}
    void SetVarImage(DataImageType *arg) {m_VarImage = arg;}

protected:
    virtual double ComputeWeightValue(unsigned int index, ImageRegionType &refPatch, ImageRegionType &movingPatch);
    virtual bool TestPatchConformity(unsigned int index, const IndexType &refIndex, const IndexType &movingIndex);

private:
    DataImagePointer m_MeanImage;
    DataImagePointer m_VarImage;

    double m_BetaParameter;
    double m_NoiseCovariance;

    double m_MeanMinThreshold;
    double m_VarMinThreshold;
};

} // end namespace anima

#include "animaNonLocalMeansPatchSearcher.hxx"
