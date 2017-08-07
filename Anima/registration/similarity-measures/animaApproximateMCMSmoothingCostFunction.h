#pragma once

#include <itkSingleValuedCostFunction.h>
#include <AnimaMCMSimilarityExport.h>
#include <animaMultiCompartmentModel.h>

namespace anima
{

class ANIMAMCMSIMILARITY_EXPORT ApproximateMCMSmoothingCostFunction : public itk::SingleValuedCostFunction
{
public:
    /** Standard class typedefs. */
    typedef ApproximateMCMSmoothingCostFunction Self;
    typedef itk::SingleValuedCostFunction Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(ApproximateMCMSmoothingCostFunction, itk::SingleValuedCostFunction)

    typedef Superclass::MeasureType MeasureType;
    typedef Superclass::DerivativeType DerivativeType;
    typedef Superclass::ParametersType ParametersType;

    typedef anima::MultiCompartmentModel MCModelType;
    typedef MCModelType::Pointer MCMPointer;
    typedef MCModelType::Vector3DType GradientType;
    typedef anima::BaseCompartment::Pointer BaseCompartmentPointer;
    typedef anima::BaseCompartment::Vector3DType Vector3DType;

    virtual MeasureType GetValue(const ParametersType &parameters) const ITK_OVERRIDE;
    virtual void GetDerivative(const ParametersType &parameters, DerivativeType &derivative) const ITK_OVERRIDE;

    void SetReferenceModels(const std::vector <MCMPointer> &refModels, const std::vector <GradientType> &gradients,
                            const std::vector <double> &bValues);
    void SetMovingModels(const std::vector<MCMPointer> &movingModels, const std::vector <GradientType> &gradients,
                         const std::vector <double> &bValues);

    void SetBValues(const std::vector <double> &val);
    void SetGradientDirections(const std::vector <GradientType> &val);

    void SetBValueWeightIndexes(const std::vector <unsigned int> &val);
    void SetSphereWeights(const std::vector <double> &val);

    itkSetMacro(ParameterScale, double)

    unsigned int GetNumberOfParameters() const ITK_OVERRIDE {return 1;}

protected:
    ApproximateMCMSmoothingCostFunction()
    {
        m_UpdatedData = false;
        m_ConstantTerm = 0;
        m_ParameterScale = 1.0e-3;
    }

    virtual ~ApproximateMCMSmoothingCostFunction() {}

private:
    ApproximateMCMSmoothingCostFunction(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    std::vector < std::vector <double> > m_ReferenceModelSignalValues;
    std::vector < std::vector <double> > m_MovingModelSignalValues;

    std::vector <unsigned int> m_BValueWeightIndexes;
    std::vector <double> m_BValues, m_SphereWeights;
    std::vector <GradientType> m_GradientDirections;

    mutable bool m_UpdatedData;
    mutable double m_ConstantTerm;
    double m_ParameterScale;
};

} // end namespace anima
