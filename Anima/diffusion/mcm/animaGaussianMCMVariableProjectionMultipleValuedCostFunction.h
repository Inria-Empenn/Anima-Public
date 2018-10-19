#pragma once

#include <itkMultipleValuedCostFunction.h>
#include <animaGaussianMCMVariableProjectionCost.h>
#include <AnimaMCMExport.h>

namespace anima
{

class ANIMAMCM_EXPORT GaussianMCMVariableProjectionMultipleValuedCostFunction : public itk::MultipleValuedCostFunction
{
public:
    /** Standard class typedefs. */
    typedef GaussianMCMVariableProjectionMultipleValuedCostFunction Self;
    typedef itk::MultipleValuedCostFunction Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(GaussianMCMVariableProjectionMultipleValuedCostFunction, itk::MultipleValuedCostFunction)

    typedef Superclass::MeasureType MeasureType;
    typedef Superclass::DerivativeType DerivativeType;
    typedef Superclass::ParametersType ParametersType;

    typedef anima::GaussianMCMVariableProjectionCost InternalCostType;
    typedef InternalCostType::Pointer InternalCostPointer;

    itkSetMacro(InternalCost, InternalCostPointer)
    itkGetConstReferenceMacro(InternalCost, InternalCostPointer)

    virtual MeasureType GetValue(const ParametersType &parameters) const ITK_OVERRIDE;
    virtual void GetDerivative(const ParametersType &parameters, DerivativeType &derivative) const ITK_OVERRIDE;

    double GetSigmaSquare();
    std::vector <double> &GetOptimalWeights();

    unsigned int GetNumberOfParameters() const ITK_OVERRIDE;
    unsigned int GetNumberOfValues() const ITK_OVERRIDE;

protected:
    GaussianMCMVariableProjectionMultipleValuedCostFunction()
    {
        m_InternalCost = ITK_NULLPTR;
    }

    virtual ~GaussianMCMVariableProjectionMultipleValuedCostFunction() {}

private:
    GaussianMCMVariableProjectionMultipleValuedCostFunction(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    mutable InternalCostPointer m_InternalCost;
};

} // end namespace anima
