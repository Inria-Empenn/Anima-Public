#pragma once
#include "animaNLMeansPatientToGroupComparisonImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkTimeProbe.h>

#include <animaNLMeansVectorPatchSearcher.h>
#include <boost/math/distributions/fisher_f.hpp>

namespace anima
{

template <class PixelScalarType>
void
NLMeansPatientToGroupComparisonImageFilter<PixelScalarType>
::BeforeThreadedGenerateData ()
{
    Superclass::BeforeThreadedGenerateData();

    // Checking consistency of the data and parameters

    unsigned int nbInputs = this->GetNumberOfIndexedInputs();
    if (nbInputs <= 0)
        itkExceptionMacro("Error: No inputs available... Exiting...");
}

template <class PixelScalarType>
void
NLMeansPatientToGroupComparisonImageFilter<PixelScalarType>
::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    typedef itk::ImageRegionConstIteratorWithIndex< InputImageType > InIteratorType;
    typedef itk::ImageRegionConstIterator < InputImageType > DatabaseIteratorType;
    typedef itk::ImageRegionIterator < OutputImageType > OutRegionIteratorType;
    typedef itk::ImageRegionConstIterator < MaskImageType > MaskRegionIteratorType;

    OutRegionIteratorType outIterator(this->GetOutput(0), outputRegionForThread);
    OutRegionIteratorType outScoreIterator(this->GetOutput(1), outputRegionForThread);
    OutRegionIteratorType outNumPatchesIterator(this->GetOutput(2), outputRegionForThread);

    MaskRegionIteratorType maskIterator (this->GetComputationMask(), outputRegionForThread);

    unsigned int numSamplesDatabase = m_DatabaseImages.size();
    std::vector <DatabaseIteratorType> databaseIterators(numSamplesDatabase);
    for (unsigned int k = 0;k < numSamplesDatabase;++k)
        databaseIterators[k] = DatabaseIteratorType(m_DatabaseImages[k],outputRegionForThread);

    InputImageType *input = const_cast<InputImageType *> (this->GetInput());
    InIteratorType patientIterator(input, outputRegionForThread);
    VectorType patientSample;

    std::vector <VectorType> databaseSamples;
    std::vector <double> databaseWeights;

    int maxAbsDisp = (int)floor((double)(m_SearchNeighborhood / m_SearchStepSize)) * m_SearchStepSize;

    typedef anima::NLMeansVectorPatchSearcher <PixelScalarType, OutputImageType> PatchSearcherType;

    PatchSearcherType patchSearcher;
    patchSearcher.SetPatchHalfSize(m_PatchHalfSize);
    patchSearcher.SetSearchStepSize(m_SearchStepSize);
    patchSearcher.SetMaxAbsDisp(maxAbsDisp);
    patchSearcher.SetInputImage(input);
    patchSearcher.SetBetaParameter(m_BetaParameter);
    patchSearcher.SetWeightThreshold(m_WeightThreshold);
    patchSearcher.SetMeanThreshold(m_MeanThreshold);
    patchSearcher.SetVarianceThreshold(m_VarianceThreshold);
    patchSearcher.SetDatabaseCovarianceDistanceAverage(m_DatabaseCovarianceDistanceAverage);
    patchSearcher.SetDatabaseCovarianceDistanceStd(m_DatabaseCovarianceDistanceStd);
    patchSearcher.SetDatabaseMeanDistanceAverage(m_DatabaseMeanDistanceAverage);
    patchSearcher.SetDatabaseMeanDistanceStd(m_DatabaseMeanDistanceStd);
    patchSearcher.SetDataMask(this->GetComputationMask());

    for (unsigned int k = 0;k < numSamplesDatabase;++k)
        patchSearcher.AddComparisonImage(m_DatabaseImages[k]);

    while (!outIterator.IsAtEnd())
    {
        if (maskIterator.Get() == 0)
        {
            outIterator.Set(0.0);
            outScoreIterator.Set(0.0);
            outNumPatchesIterator.Set(0.0);
            ++outIterator;
            ++outScoreIterator;
            ++maskIterator;
            ++patientIterator;
            ++outNumPatchesIterator;

            for (unsigned int k = 0;k < numSamplesDatabase;++k)
                ++databaseIterators[k];

            continue;
        }

        patchSearcher.UpdateAtPosition(patientIterator.GetIndex());

        databaseSamples = patchSearcher.GetDatabaseSamples();
        databaseWeights = patchSearcher.GetDatabaseWeights();

        // Add center pixels
        for (unsigned int k = 0;k < numSamplesDatabase;++k)
        {
            databaseSamples.push_back(databaseIterators[k].Get());
            databaseWeights.push_back(1.0);
        }

        double diffScore = 0;
        patientSample = patientIterator.Get();
        double pValue = 0;

        try
        {
            if (databaseSamples.size() > 1)
                pValue = this->ComputeWeightedDistanceScore(patientSample,databaseWeights,databaseSamples,diffScore);
        }
        catch (...)
        {
            itkExceptionMacro("dying in final calculations");
        }

        outIterator.Set(pValue);
        outScoreIterator.Set(diffScore);
        outNumPatchesIterator.Set(databaseSamples.size());

        this->IncrementNumberOfProcessedPoints();
        ++outIterator;
        ++outScoreIterator;
        ++outNumPatchesIterator;
        ++maskIterator;
        ++patientIterator;

        for (unsigned int k = 0;k < numSamplesDatabase;++k)
            ++databaseIterators[k];
    }
}

template <class PixelScalarType>
double
NLMeansPatientToGroupComparisonImageFilter<PixelScalarType>
::ComputeWeightedDistanceScore(itk::VariableLengthVector <double> &patientSample, std::vector < double > &databaseWeights,
                             std::vector < itk::VariableLengthVector <double> > &databaseSamples, double &diffScore)
{
    diffScore = 0;

    unsigned int nsamples = databaseSamples.size();

    unsigned int ndim = patientSample.GetSize();

    if (nsamples <= ndim)
        return 0;

    // Compute patient sample diff score
    double sumWeights = 0;
    itk::VariableLengthVector <double> databaseMean(ndim);
    databaseMean.Fill(0);

    for (unsigned int i = 0;i < nsamples;++i)
        sumWeights += databaseWeights[i];

    double sqrSumWeights = 0;
    for (unsigned int i = 0;i < nsamples;++i)
    {
        databaseWeights[i] /= sumWeights;
        sqrSumWeights += databaseWeights[i] * databaseWeights[i];
    }

    for (unsigned int i = 0;i < nsamples;++i)
    {
        databaseMean += databaseWeights[i] * databaseSamples[i];
    }

    vnl_matrix <double> covMatrixDatabase(ndim,ndim);
    covMatrixDatabase.fill(0);

    for (unsigned int i = 0;i < nsamples;++i)
    {
        for (unsigned int j = 0;j < ndim;++j)
            for (unsigned int k = j;k < ndim;++k)
                covMatrixDatabase(j,k) += databaseWeights[i]*(databaseMean[j] - databaseSamples[i][j])*(databaseMean[k] - databaseSamples[i][k]);
    }

    for (unsigned int j = 0;j < ndim;++j)
        for (unsigned int k = j;k < ndim;++k)
        {
            covMatrixDatabase(j,k) *= 1.0 / (1.0 - sqrSumWeights);
            if (j != k)
                covMatrixDatabase(k,j) = covMatrixDatabase(j,k);
        }

    vnl_matrix_inverse <double> invCovMatrixDatabase(covMatrixDatabase);
    vnl_matrix <double> covInvDatabase = invCovMatrixDatabase.inverse();

    diffScore = 0;
    for (unsigned int j = 0;j < ndim;++j)
        for (unsigned int k = j;k < ndim;++k)
        {
            double factor = 2.0;
            if (j == k)
                factor = 1;

            diffScore += factor*(patientSample[j] - databaseMean[j])*(patientSample[k] - databaseMean[k])*covInvDatabase(j,k);
        }

    double testScore = nsamples * (nsamples - ndim) * diffScore / ((nsamples *nsamples - 1) * ndim);
    boost::math::fisher_f_distribution <double> fisherDist(ndim,nsamples - ndim);
    double resVal = 1.0 - boost::math::cdf(fisherDist, testScore);

    diffScore = std::sqrt(diffScore);

    return resVal;
}

} //end namespace anima
