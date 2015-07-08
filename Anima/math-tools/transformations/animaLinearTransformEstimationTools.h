#pragma once

#include <itkAffineTransform.h>
#include <itkPoint.h>
#include <itkVector.h>
#include <itkImage.h>

namespace anima
{

// Transformation estimation utilities, many of these tools are implementations of X. Pennec PhD thesis (chap. 8)

template <class TInput, class TScalarType, unsigned int NDimensions>
        void computeLogEuclideanAverage(std::vector < vnl_matrix <TInput> > &inputTransforms, std::vector <TInput> &weights,
                                        typename itk::AffineTransform<TScalarType,NDimensions>::Pointer &resultTransform);

template <class TInput, class TScalarType, unsigned int NDimensions>
void computeTranslationLSWFromTranslations(std::vector < itk::Point<TInput,NDimensions> > &inputOrigins,
                                           std::vector < itk::Point<TInput,NDimensions> > &inputTransformed,
                                           std::vector <TInput> &weights,
                                           typename itk::AffineTransform<TScalarType,NDimensions>::Pointer &resultTransform);

template <class TInput, class TScalarType, unsigned int NDimensions>
void computeRigidLSWFromTranslations(std::vector < itk::Point<TInput,NDimensions> > &inputOrigins,
                                     std::vector < itk::Point<TInput,NDimensions> > &inputTransformed,
                                     std::vector <TInput> &weights,
                                     typename itk::AffineTransform<TScalarType,NDimensions>::Pointer &resultTransform);

template <class TInput, class TScalarType, unsigned int NDimensions>
void computeAffineLSWFromTranslations(std::vector < itk::Point<TInput,NDimensions> > &inputOrigins,
                                      std::vector < itk::Point<TInput,NDimensions> > &inputTransformed,
                                      std::vector <TInput> &weights,
                                      typename itk::AffineTransform<TScalarType,NDimensions>::Pointer &resultTransform);

// Quaternion utilities
template <class TInput, class TOutput> vnl_matrix <TOutput> computeRotationFromQuaternion(vnl_vector <TInput> eigenVector);

template <class TInput, class TOutput, unsigned int NDimensions> void pairingToQuaternion(itk::Vector<TInput,NDimensions> inputPoint,
                                                                                          itk::Vector<TInput,NDimensions> inputTransformedPoint,
                                                                                          vnl_matrix <TOutput> &outputMatrix);

}// end of namespace anima

#include "animaLinearTransformEstimationTools.hxx"
