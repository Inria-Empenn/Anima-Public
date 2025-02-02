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
itk::Point <TInput, NDimensions> computeAnisotropSimLSWFromTranslations(std::vector < itk::Point<TInput, NDimensions> > &inputOrigins,
                                                                        std::vector < itk::Point<TInput, NDimensions> > &inputTransformed,
                                                                        std::vector <TInput> &weights,
                                                                        typename itk::AffineTransform<TScalarType, NDimensions>::Pointer &resultTransform,
                                                                        vnl_matrix <double> &UMatrix);

template <class TInput, class TScalarType, unsigned int NDimensions>
itk::Point <TInput, NDimensions> computeAffineLSWFromTranslations(std::vector < itk::Point<TInput,NDimensions> > &inputOrigins,
                                                                  std::vector < itk::Point<TInput,NDimensions> > &inputTransformed,
                                                                  std::vector <TInput> &weights,
                                                                  typename itk::AffineTransform<TScalarType,NDimensions>::Pointer &resultTransform);

}// end of namespace anima

#include "animaLinearTransformEstimationTools.hxx"
