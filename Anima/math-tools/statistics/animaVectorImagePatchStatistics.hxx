#pragma once
#include "animaVectorImagePatchStatistics.h"

#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkSymmetricEigenAnalysis.h>

namespace anima
{

template <class T1, class T2, unsigned int Dimension>
unsigned int computePatchMeanAndCovariance(const itk::VectorImage <T1, Dimension> *inputImage, const itk::ImageRegion<Dimension> &patchRegion,
                                           itk::VariableLengthVector <T2> &patchMean, vnl_matrix <T2> &patchCov)
{
    unsigned int ndim = inputImage->GetNumberOfComponentsPerPixel();
    if (patchMean.GetSize() != ndim)
        patchMean.SetSize(ndim);
    patchMean.Fill(0);

    patchCov.set_size(ndim,ndim);
    patchCov.fill(0);

    typedef itk::VectorImage <T1, Dimension> VectorImageType;
    typedef itk::ImageRegionConstIteratorWithIndex< VectorImageType > InIteratorType;

    InIteratorType imageIt(inputImage,patchRegion);

    unsigned int numPixels = 0;
    while (!imageIt.IsAtEnd())
    {
        patchMean += imageIt.Get();
        ++numPixels;
        ++imageIt;
    }

    patchMean /= numPixels;
    imageIt.GoToBegin();

    itk::VariableLengthVector <T2> patchDiff(ndim);
    patchDiff.Fill(0);

    while (!imageIt.IsAtEnd())
    {
        for (unsigned int i = 0;i < ndim;++i)
            patchDiff = patchMean[i] - imageIt.Get()[i];

        for (unsigned int i = 0;i < ndim;++i)
            for (unsigned int j = i;j < ndim;++j)
            {
                double sqDiff = patchDiff[i] * patchDiff[j];
                patchCov(i,j) += sqDiff;
            }

        ++imageIt;
    }

    patchCov /= (numPixels - 1.0);
    for (unsigned int i = 0;i < ndim;++i)
        for (unsigned int j = i+1;j < ndim;++j)
            patchCov(j,i) = patchCov(i,j);

    return numPixels;
}

template <class T1, class T2, unsigned int Dimension>
void computeAverageLocalCovariance(vnl_matrix <T2> &resVariance, itk::VectorImage <T1, Dimension> *inputImage,
                                   itk::Image<unsigned char, Dimension> *maskImage, const itk::ImageRegion<Dimension> &averagingRegion,
                                   int localNeighborhood)
{
    //Generalization of what Pierrick Coupe was using in his NL-Means on scalar images

    typedef itk::VectorImage <T1, Dimension> VectorImageType;
    typedef itk::ImageRegionConstIteratorWithIndex< VectorImageType > InIteratorType;
    typedef itk::ImageRegionConstIterator < itk::Image<unsigned char, Dimension> > MaskRegionIteratorType;

    typedef itk::ImageRegion<Dimension> ImageRegionType;
    typedef typename VectorImageType::IndexType ImageIndexType;
    typedef itk::VariableLengthVector <T2> VectorType;

    ImageRegionType largestRegion = inputImage->GetLargestPossibleRegion();

    MaskRegionIteratorType maskIterator (maskImage, averagingRegion);
    InIteratorType dataIterator (inputImage, averagingRegion);

    unsigned int ndim = inputImage->GetNumberOfComponentsPerPixel();
    resVariance.set_size(ndim,ndim);
    resVariance.fill(0);

    VectorType averageLocalSignal(ndim), diffSignal(ndim);
    VectorType baseSignal(ndim);
    ImageIndexType baseIndex;

    unsigned int numEstimations = 0;

    while (!maskIterator.IsAtEnd())
    {
        if (maskIterator.Get() == 0)
        {
            ++maskIterator;
            ++dataIterator;
            continue;
        }

        baseSignal = dataIterator.Get();
        baseIndex = dataIterator.GetIndex();
        averageLocalSignal = - baseSignal;

        ImageRegionType diffRegion;
        for (unsigned int i = 0;i < Dimension;++i)
        {
            diffRegion.SetIndex(i,std::max(0,(int)baseIndex[i] - localNeighborhood));
            diffRegion.SetSize(i,std::min((unsigned int)(largestRegion.GetSize()[i] - 1),(unsigned int)(baseIndex[i] + localNeighborhood)) - diffRegion.GetIndex(i) + 1);
        }

        unsigned int numLocalPixels = diffRegion.GetSize()[0];
        for (unsigned int i = 1;i < Dimension;++i)
            numLocalPixels *= diffRegion.GetSize()[i];

        InIteratorType localIterator(inputImage,diffRegion);

        while (!localIterator.IsAtEnd())
        {
            averageLocalSignal += localIterator.Get();
            ++localIterator;
        }

        averageLocalSignal /= numLocalPixels;
        diffSignal = sqrt(numLocalPixels / (numLocalPixels + 1.0)) * (baseSignal - averageLocalSignal);

        for (unsigned int i = 0;i < ndim;++i)
            for (unsigned int j = i;j < ndim;++j)
            {
                double tmpVal = diffSignal[i] * diffSignal[j];
                resVariance(i,j) += tmpVal;
                if (i != j)
                    resVariance(j,i) += tmpVal;
            }

        ++numEstimations;
        ++maskIterator;
        ++dataIterator;
    }

    // Now divide by number of estimations and compute average variance
    resVariance /= numEstimations;
}

template <class T> double VectorCovarianceTest(vnl_matrix <T> &logRefPatchCov, vnl_matrix <T> &movingPatchCov)
{
    // Here comes a test on the distance between variances
    unsigned int ndim = logRefPatchCov.rows();

    vnl_matrix <double> eVec(ndim,ndim);
    vnl_diag_matrix <double> eVals(ndim);
    itk::SymmetricEigenAnalysis <vnl_matrix <T>, vnl_diag_matrix <double>, vnl_matrix <double> > EigenAnalysis(ndim);

    EigenAnalysis.ComputeEigenValuesAndVectors(movingPatchCov, eVals, eVec);

    for (unsigned int i = 0;i < ndim;++i)
        eVals[i] = log(eVals[i]);

    vnl_matrix <double> logMoving = eVec.transpose() * eVals * eVec;

    double varsDist = 0;
    for (unsigned int i = 0;i < ndim;++i)
        for (unsigned int j = i;j < ndim;++j)
        {
            if (i == j)
                varsDist += (logRefPatchCov(i,j) - logMoving(i,j)) * (logRefPatchCov(i,j) - logMoving(i,j));
            else
                varsDist += 2.0 * (logRefPatchCov(i,j) - logMoving(i,j)) * (logRefPatchCov(i,j) - logMoving(i,j));
        }

    varsDist = sqrt(varsDist);

    return varsDist;
}

template <class T>
double VectorMeansTest(itk::VariableLengthVector <T> &refPatchMean, itk::VariableLengthVector <T> &movingPatchMean,
                       const unsigned int &refPatchNumElts, const unsigned int &movingPatchNumElts,
                       vnl_matrix <T> &refPatchCov, vnl_matrix <T> &movingPatchCov)
{
    // Here comes a test on the distance between variances
    unsigned int ndim = refPatchMean.GetNumberOfElements();

    itk::VariableLengthVector <T> meansDiff;

    vnl_matrix <double> S_pooled(ndim,ndim), S_pooled_inv(ndim,ndim);

    S_pooled = refPatchCov * refPatchNumElts + movingPatchCov * movingPatchNumElts;
    S_pooled /= (refPatchNumElts + movingPatchNumElts - 2.0);
    S_pooled_inv = vnl_matrix_inverse <double> (S_pooled);

    meansDiff = refPatchMean - movingPatchMean;
    double distMeans = 0;

    for (unsigned int j = 0;j < ndim;++j)
        for (unsigned int k = j;k < ndim;++k)
        {
            if (j != k)
                distMeans += 2.0 * S_pooled_inv(j,k) * meansDiff[j] * meansDiff[k];
            else
                distMeans += S_pooled_inv(j,k) * meansDiff[j] * meansDiff[k];
        }

    distMeans = sqrt(distMeans);

    return distMeans;
}

} // end of namespace anima
