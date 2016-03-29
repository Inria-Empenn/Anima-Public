#pragma once

#include <vector>
#include <itkSingleValuedCostFunction.h>
#include "AnimaSHToolsExport.h"

namespace anima
{

class ANIMASHTOOLS_EXPORT ODFMaximaCostFunction :
public itk::SingleValuedCostFunction
{
public:
    /** Standard class typedefs. */
    typedef ODFMaximaCostFunction        Self;
    typedef SingleValuedCostFunction     Superclass;
    typedef itk::SmartPointer<Self>           Pointer;
    typedef itk::SmartPointer<const Self>     ConstPointer;

    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(ODFMaximaCostFunction, SingleValuedCostFunction);

    typedef Superclass::MeasureType    MeasureType;
    typedef Superclass::DerivativeType DerivativeType;
    typedef Superclass::ParametersType ParametersType;

    virtual MeasureType GetValue(const ParametersType & parameters) const ITK_OVERRIDE;
    virtual void GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const ITK_OVERRIDE;

    void SetBasisParameters(const std::vector <double> &basisPars) {m_BasisParameters = basisPars;}
    void SetODFSHOrder(unsigned int num) {m_ODFSHOrder = num;}

    virtual unsigned int GetNumberOfParameters() const ITK_OVERRIDE
    {
        return 2;
    }

protected:
    ODFMaximaCostFunction()
    {
        m_ODFSHOrder = 4;
    }

    virtual ~ODFMaximaCostFunction() {}

private:
    ODFMaximaCostFunction(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    std::vector <double> m_BasisParameters;
    unsigned int m_ODFSHOrder;
};

} // end of namespace anima


