#pragma once
#include "animaNLMeansVectorPatchSearcher.h"

#include <animaVectorImagePatchStatistics.h>
#include <itkSymmetricEigenAnalysis.h>
#include <itkImageRegionConstIterator.h>

namespace anima
{

template <class ImageScalarType, class DataImageType>
NLMeansVectorPatchSearcher <ImageScalarType, DataImageType>
::NLMeansVectorPatchSearcher()
{
    m_BetaParameter = 1.0;

    m_MeanThreshold = 0.95;
    m_VarianceThreshold = 0.5;
}

template <class ImageScalarType, class DataImageType>
void
NLMeansVectorPatchSearcher <ImageScalarType, DataImageType>
::ComputeInputProperties(const IndexType &refIndex, ImageRegionType &refPatch)
{
    unsigned int ndim = this->GetInputImage()->GetVectorLength();
    if (m_RefPatchMean.GetSize() != ndim)
    {
        m_RefPatchMean.SetSize(ndim);
        m_RefPatchCovariance.set_size(ndim,ndim);
        m_LogRefPatchCovariance.set_size(ndim,ndim);
    }

    m_RefPatchNumElements = anima::computePatchMeanAndCovariance(this->GetInputImage(),refPatch,
                                                                 m_RefPatchMean,m_RefPatchCovariance);

    CovarianceType eVec(ndim,ndim);
    vnl_diag_matrix <double> eVals(ndim);

    itk::SymmetricEigenAnalysis <vnl_matrix <double>, vnl_diag_matrix <double>, vnl_matrix <double> > EigenAnalysis(ndim);
    EigenAnalysis.ComputeEigenValuesAndVectors(m_RefPatchCovariance, eVals, eVec);

    for (unsigned int i = 0;i < ndim;++i)
        eVals[i] = std::log(eVals[i]);

    m_LogRefPatchCovariance = eVec.transpose() * eVals * eVec;

    ImageRegionType regionLocalVariance = refPatch;
    for (unsigned int i = 0;i < DataImageType::ImageDimension;++i)
    {
        regionLocalVariance.SetIndex(i,std::max(0,(int)refIndex[i] - 2));
        regionLocalVariance.SetSize(i,std::min((unsigned int)(this->GetInputImage()->GetLargestPossibleRegion().GetSize()[i] - 1),(unsigned int)(refIndex[i] + 2)) - regionLocalVariance.GetIndex(i) + 1);
    }

    anima::computeAverageLocalCovariance(m_NoiseCovariance,this->GetInputImage(),m_DataMask.GetPointer(),regionLocalVariance,2);
    m_NoiseSigma = vnl_matrix_inverse <double> (m_NoiseCovariance);
}

template <class ImageScalarType, class DataImageType>
void
NLMeansVectorPatchSearcher <ImageScalarType, DataImageType>
::ComputeComparisonProperties(unsigned int index, ImageRegionType &movingPatch)
{
    unsigned int ndim = this->GetInputImage()->GetVectorLength();
    if (m_MovingPatchMean.GetSize() != ndim)
    {
        m_MovingPatchMean.SetSize(ndim);
        m_MovingPatchCovariance.set_size(ndim,ndim);
    }

    m_MovingPatchNumElements = anima::computePatchMeanAndCovariance(this->GetComparisonImage(index),movingPatch,
                                                                    m_MovingPatchMean,m_MovingPatchCovariance);
}

template <class ImageScalarType, class DataImageType>
bool
NLMeansVectorPatchSearcher <ImageScalarType, DataImageType>
::TestPatchConformity(unsigned int index, const IndexType &refIndex, const IndexType &movingIndex)
{
    double covStdValue = m_DatabaseCovarianceDistanceStd->GetPixel(refIndex);
    double covMeanValue = m_DatabaseCovarianceDistanceAverage->GetPixel(refIndex);
    double meanStdValue = m_DatabaseMeanDistanceStd->GetPixel(refIndex);
    double meanMeanValue = m_DatabaseMeanDistanceAverage->GetPixel(refIndex);

    double varTest = anima::VectorCovarianceTest(m_LogRefPatchCovariance,m_MovingPatchCovariance);

    if (varTest > covMeanValue + m_VarianceThreshold * covStdValue)
        return false;

    double meanTestValue = anima::VectorMeansTest(m_RefPatchMean,m_MovingPatchMean,m_RefPatchNumElements,
                                                  m_MovingPatchNumElements,m_RefPatchCovariance,m_MovingPatchCovariance);

    if (meanTestValue > meanMeanValue + m_MeanThreshold * meanStdValue)
        return false;

    return true;
}

template <class ImageScalarType, class DataImageType>
double
NLMeansVectorPatchSearcher <ImageScalarType, DataImageType>
::ComputeWeightValue(unsigned int index, ImageRegionType &refPatch, ImageRegionType &movingPatch)
{
    typedef itk::ImageRegionConstIterator <RefImageType> InIteratorType;

    InIteratorType tmpIt (this->GetInputImage(), refPatch);
    InIteratorType tmpMovingIt (this->GetComparisonImage(index), movingPatch);

    VectorType tmpDiffValue;

    double weightValue = 0.0;
    unsigned int numVoxels = 0;
    unsigned int ndim = this->GetInputImage()->GetVectorLength();

    while (!tmpIt.IsAtEnd())
    {
        tmpDiffValue = tmpIt.Get() - tmpMovingIt.Get();
        for (unsigned int i = 0;i < ndim;++i)
            for (unsigned int j = i;j < ndim;++j)
            {
                if (j != i)
                    weightValue += 2.0 * m_NoiseSigma(i,j) * tmpDiffValue[i] * tmpDiffValue[j];
                else
                    weightValue += m_NoiseSigma(i,j) * tmpDiffValue[i] * tmpDiffValue[j];
            }

        ++numVoxels;
        ++tmpIt;
        ++tmpMovingIt;
    }

    weightValue = std::exp(- weightValue / (2.0 * ndim * m_BetaParameter * numVoxels));
    return weightValue;
}

} // end namespace anima
