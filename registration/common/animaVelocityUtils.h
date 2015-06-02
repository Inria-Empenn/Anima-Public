#pragma once

#include <itkStationaryVelocityFieldTransform.h>
#include <rpiDisplacementFieldTransform.h>

namespace anima
{

template <class ScalarType, unsigned int NDimensions>
void composeSVF(itk::StationaryVelocityFieldTransform <ScalarType,NDimensions> *baseTrsf,
                itk::StationaryVelocityFieldTransform <ScalarType,NDimensions> *addonTrsf);

template <class ScalarType, unsigned int NDimensions>
void GetSVFExponential(itk::StationaryVelocityFieldTransform <ScalarType,NDimensions> *baseTrsf,
                       rpi::DisplacementFieldTransform <ScalarType,NDimensions> *resultTransform,
                       bool invert);

template <class ScalarType, unsigned int NDimensions>
void composeDistortionCorrections(rpi::DisplacementFieldTransform <ScalarType,NDimensions> *baseTrsf,
                                  rpi::DisplacementFieldTransform <ScalarType,NDimensions> *positiveAddOn,
                                  rpi::DisplacementFieldTransform <ScalarType,NDimensions> *negativeAddOn,
                                  unsigned int numThreads);

} // end of namespace anima

#include "animaVelocityUtils.hxx"
