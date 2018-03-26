#include "animaGammaMixtureT2RelaxometryCostFunction.h"

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/special_functions/polygamma.hpp>

#include <vnl/vnl_matrix.h>
#include <vnl_vector_fixed.h>
#include <vnl/algo/vnl_svd.h>

namespace anima
{

GammaMixtureT2RelaxometryCostFunction::GammaMixtureT2RelaxometryCostFunction()
{
}

GammaMixtureT2RelaxometryCostFunction::MeasureType GammaMixtureT2RelaxometryCostFunction::GetValue(const ParametersType & parameters) const
{
    vnl_vector_fixed <double, 3> MeanParam;

    // First handle input parameters depending on constrained or unconstrained mode
    if (m_ConstrainedParameters)
    {
        MeanParam[0] = m_MeanParam[0];
        MeanParam[1] = parameters[0];
        MeanParam[2] = m_MeanParam[1];
    }
    else
    {
        for (unsigned int i = 0;i < 3;++i)
            MeanParam[i] = parameters[i];
    }

    m_Lambda.set_size(m_NEchoes, 3);
    m_Lambda.fill(0.0);

    /* Formulation of the LAMBDA matrix. Now on referred as 'L' for commenting */
    for (unsigned int j = 0;j < 3;++j)
    {
        double shape = MeanParam[j] * MeanParam[j] / m_VarParam[j];
        double scale = m_VarParam[j] / MeanParam[j];

        for (unsigned int k = 1;k < m_NumSteps;++k)
        {
            // The last t2 value is ignored to resolve calculation of last step issue. Shall not have worthy effect on results
            long double Gamma_Dist_Val = boost::math::gamma_p_derivative(shape, m_T2WorkingValues[k] / scale) / scale;

            for (unsigned int i = 0;i < m_NEchoes;++i)
                m_Lambda(i, j) += Gamma_Dist_Val * m_EPGSignalValues[k][i] * (m_T2WorkingValues[k] - m_T2WorkingValues[k-1]);
        }
    }

    // Get pseudo-inverse of Lambda, i.e. 'L'
    vnl_svd <double> svd(m_Lambda, 1e-08);
    m_Lambda_pinv = svd.pinverse();

    // If any weight is '0' --> Set Corresponding Lamda Col to '0'
    bool changedLambda = false;
    for (unsigned int i = 0;i < 3;++i)
    {
        double tmpWeight = 0;
        for (unsigned int j = 0;j < m_NEchoes;++j)
            tmpWeight += m_Lambda_pinv(i,j) * m_SignalValues[j];

        if (tmpWeight < 0)
        {
            changedLambda = true;
            for (unsigned int r = 0;r < m_NEchoes;++r)
                m_Lambda(r, i) = 0;
        }
    }

    //Recompute the PseudoInverse Value of Lambda
    if (changedLambda)
    {
        vnl_svd<double> svd1(m_Lambda, 1e-08);
        m_Lambda_pinv = svd1.pinverse();
    }

    // Get the Product of Lambda and its pseudo-inverse i.e. L * pinv(L)
    m_OrthoProjLambda.set_size(m_NEchoes,m_NEchoes);
    m_OrthoProjLambda.fill(0.0);

    for (unsigned int i = 0;i < m_NEchoes;++i)
    {
        m_OrthoProjLambda(i,i) = 1.0;
        for (unsigned int j = 0;j < m_NEchoes;++j)
        {
            for (unsigned int k = 0;k < 3;++k)
                m_OrthoProjLambda(i,j) -= m_Lambda(i,k) * m_Lambda_pinv(k,j);
        }
    }

    // ****** Computation of Cost Function ****** //
    m_CostFunctionVec.set_size(m_NEchoes);
    m_CostFunctionVec.fill(0);

    for (unsigned int i = 0;i < m_NEchoes;++i)
    {
        m_CostFunctionVec[i] = 0;
        for (unsigned int j = 0;j < m_NEchoes;++j)
            m_CostFunctionVec[i] += m_OrthoProjLambda(i,j) * m_SignalValues[j];
    }

    long double CostFunction_val = 0.0;
    for (unsigned int i = 0;i < m_NEchoes;++i)
        CostFunction_val += m_CostFunctionVec[i] * m_CostFunctionVec[i];

    return CostFunction_val;
}

/* main part where we calculate the derivative updates */
void GammaMixtureT2RelaxometryCostFunction::GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const
{
    unsigned int numParams = this->GetNumberOfParameters();
    derivative.SetSize(numParams);

    vnl_vector_fixed <double, 3> MeanParam;

    if (m_ConstrainedParameters)
    {
        MeanParam[0] = m_MeanParam[0];
        MeanParam[1] = parameters[0];
        MeanParam[2] = m_MeanParam[1];
    }
    else
    {
        for (unsigned int i = 0;i < 3;++i)
            MeanParam[i] = parameters[i];
    }

    /* Jacobian Computation */
    m_Jacobian_Update.set_size(m_NEchoes, numParams);
    m_Partial_Derivative.set_size(m_NEchoes, 3);
    m_DerivativeProduct.set_size(m_NEchoes, m_NEchoes);

    for (unsigned int j = 0;j < numParams;++j)
    {
        m_Partial_Derivative.fill(0.0);

        unsigned int idx = j;
        if (m_ConstrainedParameters)
            idx = 1;

        double shape = MeanParam[idx] * MeanParam[idx] / m_VarParam[idx];
        double scale = m_VarParam[idx] / MeanParam[idx];

        // Populate the matrix of differentiation with Mean-Parameter
        for (unsigned int k = 1; k < m_NumSteps; ++k)
        {
            double gammaDistValue = boost::math::gamma_p_derivative(shape, m_T2WorkingValues[k] / scale) / scale;
            double internalTerm = (2.0 * std::log(m_T2WorkingValues[k] / scale) - 2.0 * boost::math::digamma(shape) + 1.0) / scale - m_T2WorkingValues[k] / m_VarParam[idx];

            for (unsigned int r = 0; r < m_NEchoes; ++r)
            {
                double lineValue = gammaDistValue * internalTerm * m_EPGSignalValues[k][r] * (m_T2WorkingValues[k] - m_T2WorkingValues[k - 1]);

                if (!m_ConstrainedParameters)
                {
                for (unsigned int c = 0; c < 3; ++c)
                    m_Partial_Derivative(r, c) += lineValue;
                }
                else
                    m_Partial_Derivative(r, idx) += lineValue;
            }
        }

        // This is the value --> [ (orthogonal_projection)*(Partial_Derivative) ] * (PseudoInv_Lambda)
        m_DerivativeProduct = m_OrthoProjLambda * m_Partial_Derivative * m_Lambda_pinv;

        for (unsigned int i = 0; i < m_NEchoes; ++i)
        {
            m_Jacobian_Update(i,j) = 0;
            for (unsigned int k = 0;k < m_NEchoes;++k)
            {
                double tmpSum = m_DerivativeProduct(i,k) + m_DerivativeProduct(k,i);
                m_Jacobian_Update(i,j) -= tmpSum * m_SignalValues[k];
            }
        }
    }

    // *****Calculate the Derivative vector***** //
    /* Method:
        The Jacobian (say J) is --> m x 2n. where m: No.of Echoes & 2n: Number of parameters
        Cost Function Vector (say r) --> m x 1
        As per Golub (2002) p.4, the Derivative of error terms, or cost function is obtained as: 2 * (J') * r
        */
    for (unsigned int i = 0;i < numParams;++i)
    {
        derivative[i] = 0.0;
        for (unsigned int j = 0;j < m_NEchoes;++j)
            derivative[i] += 2.0 * m_Jacobian_Update(j,i) * m_CostFunctionVec[j];
    }
}

} // end namespace anima
