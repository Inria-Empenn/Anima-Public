#pragma once
#include <itkSingleValuedCostFunction.h>
#include "AnimaRelaxometryExport.h"

#include <animaNNLSOptimizer.h>
#include <vnl/vnl_matrix.h>

namespace anima
{

/**
 * \class MultiT2TikhonovCostFunction
 * \brief Cost function for estimating final solution of multi T2 estimation with Tikhonov regularization
 */
class ANIMARELAXOMETRY_EXPORT MultiT2TikhonovCostFunction :
        public itk::SingleValuedCostFunction
{
public:
    /** Standard class typedefs. */
    typedef MultiT2TikhonovCostFunction Self;
    typedef itk::SingleValuedCostFunction Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(MultiT2TikhonovCostFunction, Superclass)

    typedef Superclass::MeasureType MeasureType;
    typedef Superclass::DerivativeType DerivativeType;
    typedef Superclass::ParametersType ParametersType;
    typedef std::vector < std::complex <double> > ComplexVectorType;
    typedef std::vector <ComplexVectorType> MatrixType;

    virtual MeasureType GetValue(const ParametersType &parameters) const ITK_OVERRIDE;
    virtual void GetDerivative(const ParametersType &parameters, DerivativeType &derivative) const ITK_OVERRIDE {}

    void SetAMatrix(vnl_matrix <double> &mat) {m_AMatrix = mat;}
    void SetT2RelaxometrySignals(ParametersType &relaxoSignals) {m_T2RelaxometrySignals = relaxoSignals;}
    void SetPriorDistribution(ParametersType &prior) {m_PriorDistribution = prior;}
    itkSetMacro(ReferenceResidual, double)
    itkGetMacro(CurrentResidual, double)

    ParametersType &GetOptimizedT2Weights() {return m_OptimizedT2Weights;}
    itkGetMacro(OptimizedM0Value, double)

    unsigned int GetNumberOfParameters() const ITK_OVERRIDE
    {
        return 1;
    }

protected:
    MultiT2TikhonovCostFunction()
    {
        m_NNLSOptimizer = anima::NNLSOptimizer::New();
    }

    virtual ~MultiT2TikhonovCostFunction() {}

private:
    MultiT2TikhonovCostFunction(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    mutable anima::NNLSOptimizer::Pointer m_NNLSOptimizer;

    double m_ReferenceResidual;
    mutable ParametersType m_T2RelaxometrySignals;
    ParametersType m_PriorDistribution;
    mutable vnl_matrix <double> m_AMatrix;
    mutable ParametersType m_OptimizedT2Weights;
    mutable double m_OptimizedM0Value;
    mutable double m_CurrentResidual;
};

} // end namespace anima
