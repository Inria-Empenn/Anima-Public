#pragma once

#include <itkSingleValuedCostFunction.h>
#include <animaBaseMCMCost.h>
#include <AnimaMCMBaseExport.h>

namespace anima
{

class ANIMAMCMBASE_EXPORT MCMSingleValuedCostFunction : public itk::SingleValuedCostFunction
{
public:
    /** Standard class typedefs. */
    typedef MCMSingleValuedCostFunction Self;
    typedef itk::SingleValuedCostFunction Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(MCMSingleValuedCostFunction, itk::SingleValuedCostFunction)

    typedef Superclass::MeasureType MeasureType;
    typedef Superclass::DerivativeType DerivativeType;
    typedef Superclass::ParametersType ParametersType;

    typedef anima::BaseMCMCost InternalCostType;
    typedef InternalCostType::Pointer InternalCostPointer;

    itkSetMacro(InternalCost, InternalCostPointer)

    virtual MeasureType GetValue(const ParametersType &parameters) const ITK_OVERRIDE;
    virtual void GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const ITK_OVERRIDE;

    double GetB0Value();
    double GetSigmaSquare();

    InternalCostType::MCMPointer &GetMCMStructure() {return m_InternalCost->GetMCMStructure();}

    unsigned int GetNumberOfParameters() const ITK_OVERRIDE;

protected:
    MCMSingleValuedCostFunction()
    {
        m_InternalCost = 0;
    }

    virtual ~MCMSingleValuedCostFunction() {}

private:
    MCMSingleValuedCostFunction(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    mutable InternalCostPointer m_InternalCost;
};

} // end namespace anima
