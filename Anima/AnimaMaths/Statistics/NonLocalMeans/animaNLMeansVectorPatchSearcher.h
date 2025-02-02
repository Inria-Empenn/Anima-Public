#pragma once

#include <animaNonLocalPatchBaseSearcher.h>
#include <itkVectorImage.h>

namespace anima
{

template <class ImageScalarType, class DataImageType>
class NLMeansVectorPatchSearcher
        : public anima::NonLocalPatchBaseSearcher < itk::VectorImage <ImageScalarType,DataImageType::ImageDimension> >
{
public:
    typedef itk::VectorImage <ImageScalarType,DataImageType::ImageDimension> RefImageType;
    typedef typename RefImageType::PixelType VectorType;
    typedef vnl_matrix <double> CovarianceType;

    typedef itk::Image <unsigned char, 3> MaskImageType;
    typedef typename MaskImageType::Pointer MaskImagePointer;

    typedef typename DataImageType::Pointer DataImagePointer;
    typedef anima::NonLocalPatchBaseSearcher <RefImageType> Superclass;
    typedef typename Superclass::ImageRegionType ImageRegionType;
    typedef typename Superclass::IndexType IndexType;

    NLMeansVectorPatchSearcher();
    virtual ~NLMeansVectorPatchSearcher() {}

    void SetBetaParameter(double arg) {m_BetaParameter = arg;}
    void SetMeanThreshold(double arg) {m_MeanThreshold = arg;}
    void SetVarianceThreshold(double arg) {m_VarianceThreshold = arg;}

    void SetDataMask(MaskImageType *arg) {m_DataMask = arg;}

    void SetDatabaseCovarianceDistanceAverage(DataImagePointer &arg) {m_DatabaseCovarianceDistanceAverage = arg;}
    void SetDatabaseCovarianceDistanceStd(DataImagePointer &arg) {m_DatabaseCovarianceDistanceStd = arg;}
    void SetDatabaseMeanDistanceAverage(DataImagePointer &arg) {m_DatabaseMeanDistanceAverage = arg;}
    void SetDatabaseMeanDistanceStd(DataImagePointer &arg) {m_DatabaseMeanDistanceStd = arg;}

protected:
    virtual void ComputeInputProperties(const IndexType &refIndex, ImageRegionType &refPatch);
    virtual void ComputeComparisonProperties(unsigned int index, ImageRegionType &movingPatch);
    virtual double ComputeWeightValue(unsigned int index, ImageRegionType &refPatch, ImageRegionType &movingPatch);
    virtual bool TestPatchConformity(unsigned int index, const IndexType &refIndex, const IndexType &movingIndex);

private:
    DataImagePointer m_DatabaseCovarianceDistanceAverage;
    DataImagePointer m_DatabaseCovarianceDistanceStd;
    DataImagePointer m_DatabaseMeanDistanceAverage;
    DataImagePointer m_DatabaseMeanDistanceStd;

    MaskImagePointer m_DataMask;

    double m_BetaParameter;

    double m_MeanThreshold;
    double m_VarianceThreshold;

    // Internal variables
    VectorType m_RefPatchMean, m_MovingPatchMean;
    CovarianceType m_RefPatchCovariance, m_MovingPatchCovariance;
    unsigned int m_RefPatchNumElements, m_MovingPatchNumElements;
    CovarianceType m_LogRefPatchCovariance;
    CovarianceType m_NoiseCovariance, m_NoiseSigma;
};

} // end namespace anima

#include "animaNLMeansVectorPatchSearcher.hxx"
