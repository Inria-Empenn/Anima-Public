#pragma once
#include "animaNLMeansPatientToGroupComparisonImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkSymmetricEigenAnalysis.h>

#include <itkTimeProbe.h>
#include <animaVectorImagePatchStatistics.h>
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
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionConstIteratorWithIndex< InputImageType > InIteratorType;
    typedef itk::ImageRegionIteratorWithIndex< OutputImageType > OutRegionIteratorType;
    typedef itk::ImageRegionConstIteratorWithIndex < MaskImageType > MaskRegionIteratorType;

    OutRegionIteratorType outIterator(this->GetOutput(0), outputRegionForThread);
    OutRegionIteratorType outScoreIterator(this->GetOutput(1), outputRegionForThread);
    OutRegionIteratorType outNumPatchesIterator(this->GetOutput(2), outputRegionForThread);

    OutRegionIteratorType dbCovarianceDistanceAverageIterator(m_DatabaseCovarianceDistanceAverage,outputRegionForThread);
    OutRegionIteratorType dbCovarianceDistanceStdIterator(m_DatabaseCovarianceDistanceStd,outputRegionForThread);
    OutRegionIteratorType dbMeanDistanceAverageIterator(m_DatabaseMeanDistanceAverage,outputRegionForThread);
    OutRegionIteratorType dbMeanDistanceStdIterator(m_DatabaseMeanDistanceStd,outputRegionForThread);

    MaskRegionIteratorType maskIterator (this->GetComputationMask(), outputRegionForThread);

    unsigned int numSamplesDatabase = m_DatabaseImages.size();
    InIteratorType patientIterator(this->GetInput(0), outputRegionForThread);

    std::vector <VectorType> databaseSamples;
    std::vector <double> databaseWeights;

    unsigned int ndim = this->GetInput(0)->GetNumberOfComponentsPerPixel();
    VectorType patientVectorValue;
    VectorType tmpDiffValue;

    OutputImageRegionType largestRegionOut = this->GetOutput()->GetLargestPossibleRegion();
    OutputImageRegionType tmpBlockRegion, tmpBlockRegionMoving;
    InputImageIndexType curIndex, dispCurIndex, centerIndex;
    centerIndex.Fill(0);

    int maxAbsDisp = (int)floor((double)(m_SearchNeighborhood / m_SearchStepSize)) * m_SearchStepSize;

    CovarianceType noiseSigma, noiseCovariance;
    VectorType refPatchMean(ndim), movingPatchMean(ndim);
    CovarianceType refPatchCov(ndim,ndim,0), logRefPatchCov(ndim,ndim,0), movingPatchCov(ndim,ndim,0);
    CovarianceType eVec(ndim,ndim);
    vnl_diag_matrix <double> eVals(ndim);

    itk::SymmetricEigenAnalysis <vnl_matrix <double>, vnl_diag_matrix <double>, vnl_matrix <double> > EigenAnalysis(ndim);

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

            ++dbCovarianceDistanceAverageIterator;
            ++dbCovarianceDistanceStdIterator;
            ++dbMeanDistanceAverageIterator;
            ++dbMeanDistanceStdIterator;

            continue;
        }

        curIndex = maskIterator.GetIndex();
        databaseWeights.clear();
        databaseSamples.clear();

        for (unsigned int i = 0;i < 3;++i)
        {
            tmpBlockRegion.SetIndex(i,std::max(0,(int)curIndex[i] - (int)m_PatchHalfSize));
            tmpBlockRegion.SetSize(i,std::min((unsigned int)(largestRegionOut.GetSize()[i] - 1),(unsigned int)(curIndex[i] + m_PatchHalfSize)) - tmpBlockRegion.GetIndex(i) + 1);
        }

        OutputImageRegionType regionLocalVariance = tmpBlockRegion;
        for (unsigned int i = 0;i < 3;++i)
        {
            regionLocalVariance.SetIndex(i,std::max(0,(int)curIndex[i] - 2));
            regionLocalVariance.SetSize(i,std::min((unsigned int)(largestRegionOut.GetSize()[i] - 1),(unsigned int)(curIndex[i] + 2)) - regionLocalVariance.GetIndex(i) + 1);
        }

        InIteratorType tmpIt (this->GetInput(0),tmpBlockRegion);

        unsigned int refPatchNumElts = anima::computePatchMeanAndCovariance(this->GetInput(0),tmpBlockRegion,refPatchMean,refPatchCov);

        EigenAnalysis.ComputeEigenValuesAndVectors(refPatchCov, eVals, eVec);

        for (unsigned int i = 0;i < ndim;++i)
            eVals[i] = log(eVals[i]);

        logRefPatchCov = eVec.transpose() * eVals * eVec;

        anima::computeAverageLocalCovariance(noiseCovariance,this->GetInput(0),this->GetComputationMask(),regionLocalVariance,2);

        noiseSigma = vnl_matrix_inverse <double> (noiseCovariance);

        double covStdValue = dbCovarianceDistanceStdIterator.Get();
        double covMeanValue = dbCovarianceDistanceAverageIterator.Get();
        double meanStdValue = dbMeanDistanceStdIterator.Get();
        double meanMeanValue = dbMeanDistanceAverageIterator.Get();

        for (unsigned int k = 0;k < numSamplesDatabase;++k)
        {
            for (int dispZ = - maxAbsDisp;dispZ <= maxAbsDisp;dispZ += m_SearchStepSize)
            {
                dispCurIndex[2] = tmpBlockRegion.GetIndex()[2] + dispZ;
                int maxBlockZ = dispCurIndex[2] + tmpBlockRegion.GetSize()[2];

                if ((dispCurIndex[2] < 0)||(maxBlockZ > (int)largestRegionOut.GetSize()[2]))
                    continue;

                for (int dispY = - maxAbsDisp;dispY <= maxAbsDisp;dispY += m_SearchStepSize)
                {
                    dispCurIndex[1] = tmpBlockRegion.GetIndex()[1] + dispY;
                    int maxBlockY = dispCurIndex[1] + tmpBlockRegion.GetSize()[1];
                    if ((dispCurIndex[1] < 0)||(maxBlockY > (int)largestRegionOut.GetSize()[1]))
                        continue;

                    for (int dispX = - maxAbsDisp;dispX <= maxAbsDisp;dispX += m_SearchStepSize)
                    {
                        dispCurIndex[0] = tmpBlockRegion.GetIndex()[0] + dispX;
                        int maxBlockX = dispCurIndex[0] + tmpBlockRegion.GetSize()[0];
                        if ((dispCurIndex[0] < 0)||(maxBlockX > (int)largestRegionOut.GetSize()[0]))
                            continue;

                        tmpBlockRegionMoving.SetIndex(dispCurIndex);
                        tmpBlockRegionMoving.SetSize(tmpBlockRegion.GetSize());

                        double weightValue = 0;
                        if ((dispX != 0)||(dispY != 0)||(dispZ != 0))
                        {
                            unsigned int movingPatchNumElts = anima::computePatchMeanAndCovariance(m_DatabaseImages[k].GetPointer(),tmpBlockRegionMoving,
                                                                                                   movingPatchMean,movingPatchCov);

                            double varTest = anima::VectorCovarianceTest(logRefPatchCov,movingPatchCov);

                            if (varTest > covMeanValue + m_VarianceThreshold * covStdValue)
                                continue;

                            double meanTestValue = anima::VectorMeansTest(refPatchMean,movingPatchMean,refPatchNumElts,movingPatchNumElts,
                                                                          refPatchCov,movingPatchCov);

                            if (meanTestValue > meanMeanValue + m_MeanThreshold * meanStdValue)
                                continue;

                            InIteratorType tmpMovingIt (m_DatabaseImages[k], tmpBlockRegionMoving);
                            tmpIt.GoToBegin();

                            unsigned int numVoxels = 0;

                            while (!tmpIt.IsAtEnd())
                            {
                                tmpDiffValue = tmpIt.Get() - tmpMovingIt.Get();

                                for (unsigned int i = 0;i < ndim;++i)
                                    for (unsigned int j = i;j < ndim;++j)
                                    {
                                        if (j != i)
                                            weightValue += 2.0 * noiseSigma(i,j) * tmpDiffValue[i] * tmpDiffValue[j];
                                        else
                                            weightValue += noiseSigma(i,j) * tmpDiffValue[i] * tmpDiffValue[j];
                                    }

                                ++numVoxels;
                                ++tmpIt;
                                ++tmpMovingIt;
                            }

                            weightValue = exp(- weightValue / (2.0 * ndim * m_BetaParameter * numVoxels));
                        }
                        else
                        {
                            weightValue = 1.0;
                        }


                        if (weightValue > m_WeightThreshold)
                        {
                            databaseWeights.push_back(weightValue);

                            centerIndex = curIndex;
                            centerIndex[0] += dispX;
                            centerIndex[1] += dispY;
                            centerIndex[2] += dispZ;

                            // Getting center index value
                            databaseSamples.push_back(m_DatabaseImages[k]->GetPixel(centerIndex));
                        }
                    }
                }
            }
        }

        //std::cout << "Got all samples, total number is: " << databaseWeights.size() << std::endl;
        //exit(-1);

        double diffScore = 0;
        VectorType patientSample = patientIterator.Get();
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

        ++outIterator;
        ++outScoreIterator;
        ++outNumPatchesIterator;
        ++maskIterator;
        ++patientIterator;

        ++dbCovarianceDistanceAverageIterator;
        ++dbCovarianceDistanceStdIterator;
        ++dbMeanDistanceAverageIterator;
        ++dbMeanDistanceStdIterator;
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

    vnl_matrix <double> covInvDatabase = vnl_matrix_inverse <double> (covMatrixDatabase);

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

    diffScore = sqrt(diffScore);

    return resVal;
}

} //end namespace anima
