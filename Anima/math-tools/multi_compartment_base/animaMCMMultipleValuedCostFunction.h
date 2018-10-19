#pragma once

#include <itkMultipleValuedCostFunction.h>
#include <animaBaseMCMCost.h>
#include <AnimaMCMBaseExport.h>

namespace anima
{

class ANIMAMCMBASE_EXPORT MCMMultipleValuedCostFunction : public itk::MultipleValuedCostFunction
{
public:
    /** Standard class typedefs. */
    typedef MCMMultipleValuedCostFunction Self;
    typedef itk::MultipleValuedCostFunction Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(MCMMultipleValuedCostFunction, itk::MultipleValuedCostFunction)

    typedef Superclass::MeasureType MeasureType;
    typedef Superclass::DerivativeType DerivativeType;
    typedef Superclass::ParametersType ParametersType;

    typedef anima::BaseMCMCost InternalCostType;
    typedef InternalCostType::Pointer InternalCostPointer;

    itkSetMacro(InternalCost, InternalCostPointer)
    itkGetConstReferenceMacro(InternalCost, InternalCostPointer)

    virtual MeasureType GetValue(const ParametersType &parameters) const ITK_OVERRIDE;
    virtual void GetDerivative(const ParametersType &parameters, DerivativeType &derivative) const ITK_OVERRIDE;

    double GetSigmaSquare();

    InternalCostType::MCMPointer &GetMCMStructure() {return m_InternalCost->GetMCMStructure();}

    unsigned int GetNumberOfParameters() const ITK_OVERRIDE;
    unsigned int GetNumberOfValues() const ITK_OVERRIDE;

protected:
    MCMMultipleValuedCostFunction()
    {
        m_InternalCost = ITK_NULLPTR;
    }

    virtual ~MCMMultipleValuedCostFunction() ITK_OVERRIDE {}

private:
    MCMMultipleValuedCostFunction(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    mutable InternalCostPointer m_InternalCost;
};

} // end namespace anima
