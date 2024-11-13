#pragma once

#include <vnl/vnl_matrix.h>
#include <vector>
#include <string>

#include "AnimaMathsSHToolsExport.h"

namespace anima
{

ANIMAMATHSSHTOOLS_EXPORT void GetEulerAnglesFromRotationMatrix(vnl_matrix <double> &R, std::vector <double> &resVal);
ANIMAMATHSSHTOOLS_EXPORT void EstimateLocalODFRotationMatrix(vnl_matrix <double> &resVal, unsigned int l,
                                    double alpha, double beta, double gamma);

ANIMAMATHSSHTOOLS_EXPORT std::vector <std::vector <double> > InitializeSampleDirections(unsigned int nbTheta, unsigned int nbPhi, std::string sampleDirFileName);

ANIMAMATHSSHTOOLS_EXPORT double GetDValue(unsigned int l, int m, int mp, double angle);

} // end of namespace anima
