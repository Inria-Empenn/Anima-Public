#pragma once

#include <cmath>
#include <math.h>
#include <fstream>

#include <itkImage.h>
#include <itkVectorImage.h>

namespace anima
{

template <class ScalarType> void ComputeVoxelWeights(const typename itk::VectorImage<ScalarType, 3>::IndexType &inputIndex,
                                                     std::vector<typename itk::VectorImage<ScalarType, 3>::IndexType> &variableIndex,
                                                     std::vector<ScalarType> &wPosition, const typename itk::Image<ScalarType, 3>::SizeType &bound,
                                                     const int step, unsigned int nbVoxels);

template <class ScalarType> void RemoveNullFascicles(std::vector<ScalarType>& v, std::vector<ScalarType>& d,
                                                     std::vector<ScalarType>& kappa,
                                                     std::vector<vnl_vector<ScalarType> >& directions,
                                                     std::vector<ScalarType>& w,
                                                     std::vector<ScalarType>& wFree,
                                                     const bool optionsBoundV);

template <class ScalarType> void CreateCovarianceMatrixFromDDIParameter(const std::vector<ScalarType>& v, const std::vector<ScalarType>& d,
                                                                        const std::vector<ScalarType>& kappa, const std::vector<vnl_vector<ScalarType> >& mu,
                                                                        std::vector<vnl_matrix<ScalarType> >& Sigma);

template <class ScalarType> void DDIAveraging(std::vector<ScalarType>& v, std::vector<ScalarType>& d,
                                              std::vector<ScalarType>& kappa,
                                              std::vector<vnl_vector<ScalarType> >& directions,
                                              std::vector<ScalarType>& w,
                                              const int method, double &averageNu, double &averageDiffusivity,
                                              double &averageKappa, vnl_vector <ScalarType> &averageDirection,
                                              double &averageWeight);

template<class ScalarType> void FreeWaterDDIAveraging(std::vector < std::vector <ScalarType> > &wFree,
                                                      std::vector < std::vector <ScalarType> > &dFree,
                                                      std::vector <ScalarType> &averageDiffusivity,
                                                      std::vector <ScalarType> &averageFreeWeight);

template <class ScalarType, class VectorType>
void ComputeDistanceMatrixBetweenFascicles(std::vector<ScalarType>& v, std::vector<ScalarType>& d,
                                           std::vector<ScalarType>& kappa,
                                           std::vector<vnl_vector<ScalarType> >& directions,
                                           int method, vnl_matrix<ScalarType>& distanceMatrix);

} // end of namespace anima

#include "animaDDIAveragingTools.hxx"
