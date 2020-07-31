#pragma once

#include <vnl/vnl_matrix.h>
#include <vector>

#include <itkSingleValuedCostFunction.h>
#include "AnimaOptimizersExport.h"

namespace anima
{
/**
 * @class BLMLambdaCostFunction
 * @brief Levenberg-Marquardt lambda update cost function (phi) used
 * for bounded levenberg marquardt optimizer @BLMLambdaCostFunction
 */
class ANIMAOPTIMIZERS_EXPORT BLMLambdaCostFunction:
        public itk::SingleValuedCostFunction
{
public:
    /** Standard class typedefs. */
    typedef BLMLambdaCostFunction Self;
    typedef itk::SingleValuedCostFunction Superclass;
    typedef itk::SmartPointer <Self> Pointer;
    typedef itk::SmartPointer <const Self> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(BLMLambdaCostFunction, SingleValuedCostFunction)

    typedef Superclass::MeasureType MeasureType;
    typedef Superclass::DerivativeType DerivativeType;
    typedef Superclass::ParametersType ParametersType;

    virtual MeasureType GetValue(const ParametersType &parameters) const ITK_OVERRIDE;
    virtual void GetDerivative(const ParametersType &parameters, DerivativeType &derivative) const ITK_OVERRIDE;

    virtual unsigned int GetNumberOfParameters() const ITK_OVERRIDE {return 1;}

    itkSetMacro(JRank, unsigned int)
    itkSetMacro(DeltaParameter, double)
    itkGetMacro(SolutionInBounds, bool)

    void SetDValues(ParametersType &dVal) {m_DValues = dVal;}
    void SetInversePivotVector(std::vector <unsigned int> &invPiv) {m_InversePivotVector = invPiv;}
    void SetPivotVector(std::vector <unsigned int> &piv) {m_PivotVector = piv;}
    void SetLowerBoundsPermutted(ParametersType &lb) {m_LowerBoundsPermutted = lb;}
    void SetUpperBoundsPermutted(ParametersType &ub) {m_UpperBoundsPermutted = ub;}
    void SetPreviousParametersPermutted(ParametersType &val) {m_PreviousParametersPermutted = val;}

    void SetInputWorkMatricesAndVectorsFromQRDerivative(vnl_matrix <double> &qrDerivative,
                                                        ParametersType &qtResiduals, unsigned int rank);

    ParametersType &GetSolutionVector() {return m_SolutionVector;}

protected:
    BLMLambdaCostFunction()
    {
    }

    virtual ~BLMLambdaCostFunction() ITK_OVERRIDE {}

    bool CheckSolutionIsInBounds(ParametersType &solutionVector) const;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(BLMLambdaCostFunction);

    ParametersType m_DValues, m_InputWResiduals, m_ZeroWResiduals;
    mutable ParametersType m_PPermutted, m_PPermuttedShrunk;
    mutable ParametersType m_WResiduals;
    std::vector <unsigned int> m_InversePivotVector, m_PivotVector;
    ParametersType m_LowerBoundsPermutted, m_UpperBoundsPermutted;
    ParametersType m_PreviousParametersPermutted;
    vnl_matrix <double> m_InputWorkMatrix;
    mutable vnl_matrix <double> m_RAlphaTranspose;
    mutable vnl_matrix <double> m_WorkMatrix;
    mutable vnl_matrix <double> m_ZeroWorkMatrix;
    mutable ParametersType m_SolutionVector;
    mutable bool m_SolutionInBounds;
    double m_DeltaParameter;
    unsigned int m_JRank;
};

} // end namespace anima
