#pragma once

#include <itkSingleValuedCostFunction.h>
#include <AnimaMCMSimilarityExport.h>
#include <animaMultiCompartmentModel.h>

namespace anima
{

class ANIMAMCMSIMILARITY_EXPORT MultiTensorSmoothingCostFunction : public itk::SingleValuedCostFunction
{
public:
    /** Standard class typedefs. */
    typedef MultiTensorSmoothingCostFunction Self;
    typedef itk::SingleValuedCostFunction Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(MultiTensorSmoothingCostFunction, itk::SingleValuedCostFunction)

    typedef Superclass::MeasureType MeasureType;
    typedef Superclass::ParametersType ParametersType;
    typedef Superclass::DerivativeType DerivativeType;

    typedef anima::MultiCompartmentModel MCModelType;
    typedef MCModelType::Pointer MCMPointer;
    typedef anima::BaseCompartment::Matrix3DType TensorType;

    virtual MeasureType GetValue(const ParametersType &parameters) const ITK_OVERRIDE;
    virtual void GetDerivative(const ParametersType &parameters, DerivativeType &derivative) const ITK_OVERRIDE;

    void SetReferenceModels(const std::vector <MCMPointer> &refModels);
    void SetMovingModels(const std::vector <MCMPointer> &movingModels);

    itkSetMacro(TensorsScale, double)
    itkSetMacro(LowPassGaussianSigma, double)

    unsigned int GetNumberOfParameters() const ITK_OVERRIDE {return 1;}

protected:
    MultiTensorSmoothingCostFunction()
    {
        m_UpdatedReferenceData = false;
        m_UpdatedMovingData = false;
        m_RecomputeConstantTerm = true;

        m_ConstantTerm = 0;
        m_TensorsScale = 1000;
        m_LowPassGaussianSigma = 2000;
    }

    virtual ~MultiTensorSmoothingCostFunction() {}

    void UpdateDeterminants() const;
    void ComputeMatrixTraces(const TensorType &matrix, double &matrixTrace, double &squaredMatrixTrace) const;

private:
    MultiTensorSmoothingCostFunction(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    std::vector < std::vector <TensorType> > m_ReferenceModels;
    std::vector < std::vector <double> > m_ReferenceModelWeights;
    std::vector <double> m_ReferenceNumberOfIsotropicCompartments;
    std::vector < std::vector <TensorType> > m_MovingModels;
    std::vector < std::vector <double> > m_MovingModelWeights;

    mutable std::vector < std::vector <double> > m_ReferenceReferenceDeterminants, m_ReferenceMovingDeterminants, m_MovingMovingDeterminants;
    mutable std::vector < std::vector <double> > m_ReferenceReferenceTraces, m_ReferenceMovingTraces;
    mutable std::vector < std::vector <double> > m_ReferenceReferenceTraceDifferences, m_ReferenceMovingTraceDifferences;

    mutable bool m_UpdatedReferenceData;
    mutable bool m_UpdatedMovingData;
    mutable bool m_RecomputeConstantTerm;

    mutable double m_ConstantTerm;
    double m_TensorsScale;
    double m_LowPassGaussianSigma;
};

} // end namespace anima
