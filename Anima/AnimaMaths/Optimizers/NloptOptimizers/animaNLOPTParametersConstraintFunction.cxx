#include <animaNLOPTParametersConstraintFunction.h>

namespace anima
{

double
NLOPTParametersConstraintFunction
::GetConstraintValue(unsigned int n, const double *x, double *grad, void *data)
{
    ConstraintDataType *dataCast = (ConstraintDataType *)data;

    return dataCast->constraintPointer->InternalComputeConstraint(n, x, grad);
}

} // end namespace anima
