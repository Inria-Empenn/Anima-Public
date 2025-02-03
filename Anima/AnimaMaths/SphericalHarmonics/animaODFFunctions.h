#pragma once

#include <vnl/vnl_matrix.h>
#include <vector>
#include <string>

#include <libAnimaMathsExport.h>

namespace anima
{

LIBANIMAMATHS_EXPORT void GetEulerAnglesFromRotationMatrix(vnl_matrix <double> &R, std::vector <double> &resVal);
LIBANIMAMATHS_EXPORT void EstimateLocalODFRotationMatrix(vnl_matrix <double> &resVal, unsigned int l,
                                                         double alpha, double beta, double gamma);

LIBANIMAMATHS_EXPORT std::vector <std::vector <double> > InitializeSampleDirections(unsigned int nbTheta, unsigned int nbPhi, std::string sampleDirFileName);

LIBANIMAMATHS_EXPORT double GetDValue(unsigned int l, int m, int mp, double angle);

} // end of namespace anima
