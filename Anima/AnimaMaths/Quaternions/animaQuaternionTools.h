#pragma once

namespace anima
{

template <class TInput, class TOutput> vnl_matrix <TOutput> computeRotationFromQuaternion(vnl_vector <TInput> eigenVector);

template <class TInput, class TOutput, unsigned int NDimensions>
void pairingToQuaternion(const vnl_vector_fixed <TInput,NDimensions> &inputPoint, const vnl_vector_fixed <TInput,NDimensions> &inputTransformedPoint,
                         vnl_matrix <TOutput> &outputMatrix);

template <class TInput, class TOutput, unsigned int NDimensions>
void pairingToQuaternion(const itk::Vector <TInput,NDimensions> &inputPoint, const itk::Vector <TInput,NDimensions> &inputTransformedPoint,
                         vnl_matrix <TOutput> &outputMatrix);

template <class PointType, class TOutput>
void pairingToQuaternion(const PointType &inputPoint, const PointType &inputTransformedPoint, vnl_matrix <TOutput> &outputMatrix, unsigned int ndim);

template <class TInput, class TOutput, unsigned int NDimensions>
void pairingToQuaternionScalingsDerivative(const vnl_vector_fixed <TInput, NDimensions> &inputPoint, const vnl_vector_fixed <TInput, NDimensions> &inputTransformedPoint,
                                           vnl_matrix <TOutput> &outputMatrix, const int &dimScal);

template <class TInput, class TOutput, unsigned int NDimensions>
void pairingToQuaternionScalingsDerivative(const itk::Vector <TInput, NDimensions> &inputPoint, const itk::Vector <TInput, NDimensions> &inputTransformedPoint,
                                           vnl_matrix <TOutput> &outputMatrix, const int &dimScal);

template <class PointType, class TOutput>
void pairingToQuaternionScalingsDerivative(const PointType &inputPoint, const PointType &inputTransformedPoint, vnl_matrix <TOutput> &outputMatrix, unsigned int ndim, const int &dimScal);

}// end of namespace anima

#include "animaQuaternionTools.hxx"
