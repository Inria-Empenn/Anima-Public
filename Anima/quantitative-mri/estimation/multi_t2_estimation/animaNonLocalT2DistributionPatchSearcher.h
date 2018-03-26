#pragma once

#include <animaNonLocalPatchBaseSearcher.h>

namespace anima
{

template <class ImageType, class DataImageType>
class NonLocalT2DistributionPatchSearcher : public anima::NonLocalPatchBaseSearcher <ImageType>
{
public:
    typedef typename DataImageType::Pointer DataImagePointer;
    typedef typename DataImageType::ConstPointer DataImageConstPointer;
    typedef anima::NonLocalPatchBaseSearcher <ImageType> Superclass;
    typedef typename Superclass::ImageRegionType ImageRegionType;
    typedef typename Superclass::IndexType IndexType;

    NonLocalT2DistributionPatchSearcher();
    virtual ~NonLocalT2DistributionPatchSearcher() {}

    void SetBetaParameter(double arg) {m_BetaParameter = arg;}
    void SetNoiseCovariances(std::vector <double> &arg) {m_NoiseCovariances = arg;}
    void SetMeanMinThreshold(double arg) {m_MeanMinThreshold = arg;}
    void SetVarMinThreshold(double arg) {m_VarMinThreshold = arg;}

    void AddPatchTestImage(const DataImageType *arg) {m_PatchTestImages.push_back(arg);}
    void AddMeanImage(DataImageType *arg) {m_MeanImages.push_back(arg);}
    void AddVarImage(DataImageType *arg) {m_VarImages.push_back(arg);}

protected:
    virtual double ComputeWeightValue(unsigned int index, ImageRegionType &refPatch, ImageRegionType &movingPatch);
    virtual bool TestPatchConformity(unsigned int index, const IndexType &refIndex, const IndexType &movingIndex);

private:
    std::vector <DataImageConstPointer> m_PatchTestImages;
    std::vector <DataImagePointer> m_MeanImages;
    std::vector <DataImagePointer> m_VarImages;

    double m_BetaParameter;
    std::vector <double> m_NoiseCovariances;

    double m_MeanMinThreshold;
    double m_VarMinThreshold;
};

} // end namespace anima

#include "animaNonLocalT2DistributionPatchSearcher.hxx"
