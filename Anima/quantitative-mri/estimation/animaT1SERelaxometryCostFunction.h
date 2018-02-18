#pragma once

#include <itkSingleValuedCostFunction.h>
#include "AnimaRelaxometryExport.h"

namespace anima
{
    
class ANIMARELAXOMETRY_EXPORT T1SERelaxometryCostFunction :
public itk::SingleValuedCostFunction
{
public:
    /** Standard class typedefs. */
    typedef T1SERelaxometryCostFunction Self;
    typedef itk::SingleValuedCostFunction   Superclass;
    typedef itk::SmartPointer<Self>         Pointer;
    typedef itk::SmartPointer<const Self>   ConstPointer;

    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(T1SERelaxometryCostFunction, Superclass)

    typedef Superclass::MeasureType MeasureType;
    typedef Superclass::DerivativeType DerivativeType;
    typedef Superclass::ParametersType ParametersType;
    typedef std::vector < std::complex <double> > ComplexVectorType;
    typedef std::vector <ComplexVectorType> MatrixType;

    virtual MeasureType GetValue(const ParametersType & parameters) const ITK_OVERRIDE;
    virtual void GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const ITK_OVERRIDE;

    void SetRelaxometrySignals(std::vector <double> & relaxoSignals) {m_RelaxometrySignals = relaxoSignals;}
    void SetTRValues(std::vector <double> & trValues) {m_TRValues = trValues;}

    unsigned int GetNumberOfParameters() const ITK_OVERRIDE
    {
        return 2;
    }

protected:
    T1SERelaxometryCostFunction()
    {
    }

    virtual ~T1SERelaxometryCostFunction() {}

private:
    T1SERelaxometryCostFunction(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    std::vector <double> m_RelaxometrySignals;
    std::vector <double> m_TRValues;
};

} // end namespace anima
