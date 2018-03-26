#include <vector>
#include <itkSingleValuedCostFunction.h>
#include <vnl/vnl_matrix.h>

#include "AnimaRelaxometryExport.h"

namespace anima
{
class ANIMARELAXOMETRY_EXPORT GammaMixtureT2RelaxometryCostFunction :
        public itk::SingleValuedCostFunction
{
public:
    /** Standard class typedefs. */
    typedef GammaMixtureT2RelaxometryCostFunction Self;
    typedef SingleValuedCostFunction Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(GammaMixtureT2RelaxometryCostFunction, SingleValuedCostFunction)

    typedef Superclass::MeasureType    MeasureType;
    typedef Superclass::DerivativeType DerivativeType;
    typedef Superclass::ParametersType ParametersType;

    /**
     * The measure type shall be used for computing the cost function value to observe convergence
     * Parameters are set as {theta_1,theta_2,theta_3} for unconstrained estimation
     * Or {theta_2} for constrained
     */
    virtual MeasureType GetValue(const ParametersType & parameters) const ITK_OVERRIDE;
    virtual void GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const ITK_OVERRIDE;

    itkSetMacro(NEchoes,int)
    itkSetMacro(NumSteps, unsigned int)

    void SetMeanParam(std::vector <double> &input) {m_MeanParam = input;}
    void SetVarParam(std::vector <double> &input) {m_VarParam = input;}

    void SetT2WorkingValues(std::vector <double> &input) {m_T2WorkingValues = input;}
    void SetEPGSignalValues(std::vector < std::vector <double> > &input) {m_EPGSignalValues = input;}
    void SetSignalValues(std::vector <double> &input) {m_SignalValues = input;}

    unsigned int GetNumberOfParameters() const ITK_OVERRIDE
    {
        if (m_ConstrainedParameters)
            return 1;
        else
            return 3;
    }

    itkSetMacro(ConstrainedParameters, bool)

    vnl_matrix <double> &GetLambda_pinv() {return m_Lambda_pinv;}

protected:
    GammaMixtureT2RelaxometryCostFunction();
    virtual ~GammaMixtureT2RelaxometryCostFunction() {}

private:
    unsigned int m_NEchoes;

    bool m_ConstrainedParameters;

    std::vector < std::vector <double> > m_EPGSignalValues;
    std::vector<double> m_T2WorkingValues;
    unsigned int m_NumSteps;
    std::vector<double> m_SignalValues;
    std::vector<double> m_VarParam;

    mutable std::vector <double> m_MeanParam;
    mutable vnl_matrix <double> m_OrthoProjLambda; // --> (I - L * pinv(L)) {where I-->Identity Matrix}, i.e. Orthogonal Projection of Lambda
    mutable vnl_matrix <double> m_Lambda;
    mutable vnl_matrix <double> m_Lambda_pinv;
    mutable vnl_vector <double> m_CostFunctionVec;

    // Internal work variables for derivative
    mutable vnl_matrix <double> m_Partial_Derivative, m_Jacobian_Update;
    mutable vnl_matrix <double> m_DerivativeProduct;
};

}
