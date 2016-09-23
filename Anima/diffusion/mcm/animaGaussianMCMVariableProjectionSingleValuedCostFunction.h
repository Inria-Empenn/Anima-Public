#pragma once

#include <itkSingleValuedCostFunction.h>
#include <animaGaussianMCMVariableProjectionCost.h>
#include <AnimaMCMExport.h>

namespace anima
{

class ANIMAMCM_EXPORT GaussianMCMVariableProjectionSingleValuedCostFunction : public itk::SingleValuedCostFunction
{
public:
    /** Standard class typedefs. */
    typedef GaussianMCMVariableProjectionSingleValuedCostFunction Self;
    typedef itk::SingleValuedCostFunction Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(GaussianMCMVariableProjectionSingleValuedCostFunction, itk::SingleValuedCostFunction)

    typedef Superclass::MeasureType MeasureType;
    typedef Superclass::DerivativeType DerivativeType;
    typedef Superclass::ParametersType ParametersType;

    typedef anima::GaussianMCMVariableProjectionCost InternalCostType;
    typedef InternalCostType::Pointer InternalCostPointer;

    itkSetMacro(InternalCost, InternalCostPointer)

    virtual MeasureType GetValue(const ParametersType &parameters) const ITK_OVERRIDE;
    virtual void GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const ITK_OVERRIDE;

    double GetB0Value();
    double GetSigmaSquare();
    std::vector <double> &GetOptimalWeights();

    unsigned int GetNumberOfParameters() const ITK_OVERRIDE;

protected:
    GaussianMCMVariableProjectionSingleValuedCostFunction()
    {
        m_InternalCost = 0;
    }

    virtual ~GaussianMCMVariableProjectionSingleValuedCostFunction() {}

private:
    GaussianMCMVariableProjectionSingleValuedCostFunction(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    mutable InternalCostPointer m_InternalCost;
};

} // end namespace anima
