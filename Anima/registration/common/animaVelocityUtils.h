#pragma once

#include <itkStationaryVelocityFieldTransform.h>
#include <rpiDisplacementFieldTransform.h>
#include <itkMultiThreader.h>

namespace anima
{

template <class ScalarType, unsigned int NDimensions>
void composeSVF(itk::StationaryVelocityFieldTransform <ScalarType,NDimensions> *baseTrsf,
                itk::StationaryVelocityFieldTransform <ScalarType,NDimensions> *addonTrsf,
                unsigned int numThreads = itk::MultiThreader::GetGlobalDefaultNumberOfThreads());

template <class ScalarType, unsigned int NDimensions>
void GetSVFExponential(itk::StationaryVelocityFieldTransform <ScalarType,NDimensions> *baseTrsf,
                       rpi::DisplacementFieldTransform <ScalarType,NDimensions> *resultTransform,
                       bool invert);

/**
 * Compose distortion correction opposite updates, ensures opposite symmetry
 * baseTrsf is replaced by the result !
 */
template <class ScalarType, unsigned int NDimensions>
void composeDistortionCorrections(typename rpi::DisplacementFieldTransform <ScalarType,NDimensions>::Pointer &baseTrsf,
                                  typename rpi::DisplacementFieldTransform <ScalarType,NDimensions>::Pointer &positiveAddOn,
                                  typename rpi::DisplacementFieldTransform <ScalarType,NDimensions>::Pointer &negativeAddOn,
                                  unsigned int numThreads =  itk::MultiThreader::GetGlobalDefaultNumberOfThreads());

} // end of namespace anima

#include "animaVelocityUtils.hxx"
