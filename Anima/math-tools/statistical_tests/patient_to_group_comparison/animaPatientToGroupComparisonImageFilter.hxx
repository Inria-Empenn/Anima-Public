#pragma once
#include "animaPatientToGroupComparisonImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkSymmetricEigenAnalysis.h>

#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/chi_squared.hpp>

#include <itkTimeProbe.h>

namespace anima
{

template <class PixelScalarType>
void
PatientToGroupComparisonImageFilter<PixelScalarType>
::BeforeThreadedGenerateData ()
{
    Superclass::BeforeThreadedGenerateData();
    // Checking consistency of the data and parameters

    unsigned int nbInputs = this->GetNumberOfIndexedInputs();

    if (nbInputs <= 0)
        itkExceptionMacro("Error: No inputs available...");

    unsigned int ndim = this->GetInput(0)->GetNumberOfComponentsPerPixel();

    if (m_NumEigenValuesPCA > ndim)
        m_NumEigenValuesPCA = ndim;

    if (m_DatabaseImages.size() <= this->GetInput(0)->GetNumberOfComponentsPerPixel())
        itkExceptionMacro("Error: Not enough inputs available...");
}

template <class PixelScalarType>
void
PatientToGroupComparisonImageFilter<PixelScalarType>
::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    typedef itk::ImageRegionConstIteratorWithIndex< InputImageType > InIteratorType;

    typedef itk::ImageRegionIteratorWithIndex< OutputImageType > OutRegionIteratorType;
    typedef itk::ImageRegionConstIteratorWithIndex < MaskImageType > MaskRegionIteratorType;

    OutRegionIteratorType outIterator(this->GetOutput(0), outputRegionForThread);
    OutRegionIteratorType outPValIterator(this->GetOutput(1), outputRegionForThread);
    MaskRegionIteratorType maskIterator (this->GetComputationMask(), outputRegionForThread);

    unsigned int numSamplesDatabase = m_DatabaseImages.size();
    std::vector < InIteratorType > databaseIterators;
    InIteratorType patientIterator(this->GetInput(0), outputRegionForThread);

    for (unsigned int i = 0;i < numSamplesDatabase;++i)
        databaseIterators.push_back(InIteratorType(m_DatabaseImages[i],outputRegionForThread));

    unsigned int ndim = this->GetInput(0)->GetNumberOfComponentsPerPixel();
    CovMatrixType tmpCovMatrix(ndim,ndim), tmpCovMatrixInv;
    VectorType tmpMean(ndim);

    std::vector < VectorType > databaseValues(numSamplesDatabase);
    VectorType patientVectorValue;

    while (!outIterator.IsAtEnd())
    {
        patientVectorValue = patientIterator.Get();
        ndim = patientVectorValue.GetSize();

        if ((maskIterator.Get() == 0)||(this->isZero(patientVectorValue)))
        {
            outIterator.Set(0.0);
            outPValIterator.Set(0.0);
            ++outIterator;
            ++outPValIterator;
            ++maskIterator;
            ++patientIterator;

            for (unsigned int i = 0;i < numSamplesDatabase;++i)
                ++databaseIterators[i];

            continue;
        }

        databaseValues.resize(numSamplesDatabase);

        unsigned int numEffectiveSamples = 0;
        for (unsigned int i = 0;i < numSamplesDatabase;++i)
        {
            databaseValues[numEffectiveSamples] = databaseIterators[i].Get();
            if (!this->isZero(databaseValues[numEffectiveSamples]))
                numEffectiveSamples++;
        }

        databaseValues.resize(numEffectiveSamples);

        ndim = this->SampleFromDiffusionModels(databaseValues,patientVectorValue);
        unsigned int ndim_afterpca = ndim;

        if ((m_NumEigenValuesPCA < ndim)||(m_ExplainedRatio < 1))
            ndim_afterpca = this->GetPCAVectorsFromData(databaseValues,patientVectorValue);

        if (numEffectiveSamples <= ndim_afterpca)
        {
            outIterator.Set(0.0);
            outPValIterator.Set(0.0);
            ++outIterator;
            ++outPValIterator;
            ++maskIterator;
            ++patientIterator;

            for (unsigned int i = 0;i < numSamplesDatabase;++i)
                ++databaseIterators[i];

            continue;
        }

        tmpCovMatrix.set_size(ndim_afterpca,ndim_afterpca);
        tmpCovMatrix.fill(0);
        tmpMean.SetSize(ndim_afterpca);
        tmpMean.Fill(0);

        for (unsigned int i = 0;i < numEffectiveSamples;++i)
        {
            for (unsigned int j = 0;j < ndim_afterpca;++j)
                tmpMean[j] += databaseValues[i][j];
        }

        tmpMean /= numEffectiveSamples;

        for (unsigned int i = 0;i < numEffectiveSamples;++i)
        {
            for (unsigned int j = 0;j < ndim_afterpca;++j)
                for (unsigned int k = j;k < ndim_afterpca;++k)
                    tmpCovMatrix(j,k) += (databaseValues[i][j] - tmpMean[j])*(databaseValues[i][k] - tmpMean[k]);
        }

        for (unsigned int j = 0;j < ndim_afterpca;++j)
            for (unsigned int k = j;k < ndim_afterpca;++k)
            {
                tmpCovMatrix(j,k) /= (numEffectiveSamples - 1.0);
                if (j != k)
                    tmpCovMatrix(k,j) = tmpCovMatrix(j,k);
            }

        tmpCovMatrixInv = vnl_matrix_inverse <double> (tmpCovMatrix);

        double resValue = 0;

        for (unsigned int i = 0;i < ndim_afterpca;++i)
            for (unsigned int j = i;j < ndim_afterpca;++j)
            {
                // This is because we are working only on triangular superior for being faster
                // The input data has to be in the right log-vector form
                double factor = 2.0;
                if (i == j)
                    factor = 1;

                resValue += factor*tmpCovMatrixInv(i,j)*(patientVectorValue[i] - tmpMean[i])*(patientVectorValue[j] - tmpMean[j]);
            }

        switch (m_StatisticalTestType)
        {
            case FISHER:
            {
                double testScore = numEffectiveSamples * (numEffectiveSamples - ndim_afterpca) * resValue / ((numEffectiveSamples * numEffectiveSamples - 1) * ndim_afterpca);
                boost::math::fisher_f_distribution <double> fisherDist(ndim_afterpca,numEffectiveSamples - ndim_afterpca);
                outPValIterator.Set (1.0 - boost::math::cdf(fisherDist, testScore));
                break;
            }

            case CHI_SQUARE:
            default:
            {
                boost::math::chi_squared_distribution <double> chiDist(ndim_afterpca);
                outPValIterator.Set (1.0 - boost::math::cdf(chiDist, resValue));
                break;
            }
        }

        resValue = sqrt(resValue);
        // If scalar values, put a sign on out z-score
        if (ndim == 1)
        {
            if (tmpMean[0] > patientVectorValue[0])
                resValue *= -1;
        }

        outIterator.Set(resValue);
        ++outIterator;
        ++outPValIterator;
        ++maskIterator;
        ++patientIterator;

        for (unsigned int i = 0;i < numSamplesDatabase;++i)
            ++databaseIterators[i];
    }
}


template <class PixelScalarType>
unsigned int
PatientToGroupComparisonImageFilter<PixelScalarType>
::GetPCAVectorsFromData(std::vector < itk::VariableLengthVector <double> > &databaseVectors,
                                                                        itk::VariableLengthVector <double> &patientVectorValue)
{
    unsigned int refNDim = patientVectorValue.GetSize();
    vnl_matrix <double> allCovarianceMatrix(refNDim,refNDim);

    for (unsigned int i = 0;i < refNDim;++i)
        for (unsigned int j = 0;j < refNDim;++j)
            allCovarianceMatrix(i,j) = 0;

    unsigned int numData = databaseVectors.size();
    std::vector <double> dataMean(refNDim,0);

    for (unsigned int i = 0;i < numData;++i)
        for (unsigned int j = 0;j < refNDim;++j)
            dataMean[j] += databaseVectors[i][j];

    for (unsigned int j = 0;j < refNDim;++j)
        dataMean[j] /= numData;

    for (unsigned int i = 0;i < numData;++i)
        for (unsigned int j = 0;j < refNDim;++j)
            for (unsigned int k = j;k < refNDim;++k)
                allCovarianceMatrix(j,k) += (databaseVectors[i][j] - dataMean[j])*(databaseVectors[i][k] - dataMean[k]);

    for (unsigned int j = 0;j < refNDim;++j)
        for (unsigned int k = j;k < refNDim;++k)
        {
            allCovarianceMatrix(j,k) /= (numData - 1.0);
            if (j != k)
                allCovarianceMatrix(k,j) = allCovarianceMatrix(j,k);
        }

    vnl_matrix <double> eigenVectors(refNDim,refNDim);
    vnl_diag_matrix <double> eigenValues(refNDim);

    itk::SymmetricEigenAnalysis <vnl_matrix <double>, vnl_diag_matrix <double>, vnl_matrix <double> > EigenAnalysis(refNDim);
    EigenAnalysis.SetOrderEigenValues(true);
    EigenAnalysis.ComputeEigenValuesAndVectors(allCovarianceMatrix,eigenValues,eigenVectors);

    double sumEigVal = 0;
    for (unsigned int i = 0;i < refNDim;++i)
        sumEigVal += eigenValues[i];

    unsigned int outNDim = 0;
    double cumulEigs = 0;
    while (cumulEigs/sumEigVal < m_ExplainedRatio)
    {
        cumulEigs += eigenValues[refNDim - 1 - outNDim];
        ++outNDim;
    }

    if (outNDim < m_NumEigenValuesPCA)
        outNDim = m_NumEigenValuesPCA;

    vnl_matrix <double> basisMatrix(outNDim,refNDim);
    for (unsigned int i = 0;i < outNDim;++i)
        for (unsigned int j = 0;j < refNDim;++j)
            basisMatrix(i,j) = eigenVectors(refNDim - 1 - i,j);

    itk::VariableLengthVector <double> tmpVec, resData(outNDim);
    for (unsigned int i = 0;i < numData;++i)
    {
        tmpVec = databaseVectors[i];
        for (unsigned int j = 0;j < outNDim;++j)
        {
            resData[j] = 0;
            for (unsigned int k = 0;k < refNDim;++k)
                resData[j] += basisMatrix(j,k)*(tmpVec[k] - dataMean[k]);
        }

        databaseVectors[i] = resData;
    }

    for (unsigned int j = 0;j < outNDim;++j)
    {
        resData[j] = 0;
        for (unsigned int k = 0;k < refNDim;++k)
            resData[j] += basisMatrix(j,k)*(patientVectorValue[k] - dataMean[k]);
    }

    patientVectorValue = resData;

    return outNDim;
}

template <class PixelScalarType>
bool
PatientToGroupComparisonImageFilter<PixelScalarType>
::isZero(const VectorType &vec)
{
    unsigned int ndim = vec.Size();

    for (unsigned int i = 0;i < ndim;++i)
    {
        if (vec[i] != 0)
            return false;
    }

    return true;
}

} //end namespace anima
