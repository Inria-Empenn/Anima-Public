#pragma once

#include <itkObject.h>
#include <itkSmartPointer.h>

#include <AnimaOptimizersExport.h>

namespace anima
{

class ANIMAOPTIMIZERS_EXPORT NLOPTParametersConstraintFunction : public itk::Object
{
public:
    /** Standard class typedefs. */
    typedef NLOPTParametersConstraintFunction Self;
    typedef itk::Object Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkTypeMacro(NLOPTParametersConstraintFunction, itk::Object)

    typedef struct
    {
        Self *constraintPointer;
    } ConstraintDataType;

    static double GetConstraintValue(unsigned int n, const double *x, double *grad, void *data);

    itkSetMacro(Tolerance, double)
    itkGetMacro(Tolerance, double)

    itkGetMacro(AdditionalData, ConstraintDataType *)

protected:
    NLOPTParametersConstraintFunction()
    {
        m_Tolerance = 1.0e-8;
        m_AdditionalData = new ConstraintDataType;
        m_AdditionalData->constraintPointer = this;
    }

    virtual ~NLOPTParametersConstraintFunction()
    {
        delete m_AdditionalData;
    }

    virtual double InternalComputeConstraint(unsigned int numParameters, const double *dataValue, double *gradValue) = 0;

private:
    NLOPTParametersConstraintFunction(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    double m_Tolerance;
    ConstraintDataType *m_AdditionalData;
};

} // end namespace anima
