#pragma once

#include <itkVectorImage.h>

namespace anima
{

//! Computes the average and covariance matrix from a patch of a vector image
template <class T1, class T2, unsigned int Dimension>
unsigned int computePatchMeanAndCovariance(const itk::VectorImage <T1, Dimension> *inputImage, const itk::ImageRegion<Dimension> &patchRegion,
                                           itk::VariableLengthVector <T2> &patchMean, vnl_matrix <T2> &patchCov);

//! Noise estimation for a patch of a vector image
template <class T1, class T2, unsigned int Dimension>
void computeAverageLocalCovariance(vnl_matrix <T2> &resVariance, itk::VectorImage <T1, Dimension> *inputImage,
                                   itk::Image<unsigned char, Dimension> *maskImage, const itk::ImageRegion<Dimension> &averagingRegion,
                                   int localNeighborhood);

//! Test if covariance matrices are different (returns distance)
template <class T> double VectorCovarianceTest(vnl_matrix <T> &logRefPatchCov, vnl_matrix <T> &movingPatchCov);

//! Test if vector means are different (returns distance)
template <class T> double VectorMeansTest(itk::VariableLengthVector <T> &refPatchMean, itk::VariableLengthVector <T> &movingPatchMean,
                                          const unsigned int &refPatchNumElts, const unsigned int &movingPatchNumElts,
                                          vnl_matrix <T> &refPatchCov, vnl_matrix <T> &movingPatchCov);

} // end of namespace anima

#include "animaVectorImagePatchStatistics.hxx"
