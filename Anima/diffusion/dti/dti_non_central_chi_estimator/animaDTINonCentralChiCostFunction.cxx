#include "animaDTINonCentralChiCostFunction.h"
#include <animaBesselFunctions.h>

namespace anima
{
    DTINonCentralChiCostFunction::DTINonCentralChiCostFunction()
    {
        m_RawSignal.clear();
        m_NumberOfParameters = 7;
        m_NumberOfCoils = 1;

        m_UseB0Value = false;
        m_B0Value = 1;
    }

    DTINonCentralChiCostFunction::MeasureType DTINonCentralChiCostFunction::GetValue(const ParametersType & parameters) const
    {
        unsigned int nbRawSignals = m_RawSignal.size();
        double costValue = - log(m_Sigma);
        costValue *= nbRawSignals;

        for (unsigned int i = 0;i < nbRawSignals;++i)
        {
            double tmpVal = 0;

            double y = m_RawSignal[i];
            tmpVal += m_NumberOfCoils * log(y);

            double modelValue = 0;
            for (unsigned int j = 0;j < m_NumberOfParameters;++j)
                modelValue += m_DesignMatrix(i,j) * parameters[j];

            modelValue = exp(modelValue);
            if (m_UseB0Value)
                modelValue *= m_B0Value;

            tmpVal -= (m_NumberOfCoils - 1) * log(modelValue);
            tmpVal -= (y * y + modelValue * modelValue) / (2.0 * m_Sigma);

            y *= modelValue / m_Sigma;
            tmpVal += anima::log_bessel_i(m_NumberOfCoils - 1,y);

            costValue += tmpVal;
        }

        return costValue;
    }

    void DTINonCentralChiCostFunction::GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const
    {
        derivative = DerivativeType( m_NumberOfParameters );
        derivative.Fill(0);

        std::vector <double> nuVector(m_RawSignal.size());
        vnl_matrix <double> nuDerivativeMatrix(m_RawSignal.size(),m_NumberOfParameters);

        for (unsigned int i = 0;i < m_RawSignal.size();++i)
        {
            double modelValue = 0;
            for (unsigned int j = 0;j < m_NumberOfParameters;++j)
                modelValue += m_DesignMatrix(i,j) * parameters[j];

            modelValue = exp(modelValue);
            if (m_UseB0Value)
                modelValue *= m_B0Value;

            nuVector[i] = modelValue;

            for (unsigned int j = 0;j < m_NumberOfParameters;++j)
                nuDerivativeMatrix(i,j) = m_DesignMatrix(i,j) * modelValue;
        }

        for (unsigned int i = 0;i < m_RawSignal.size();++i)
        {
            double dnuValue = - nuVector[i];
            double y = m_RawSignal[i];

            double insideValue = y * nuVector[i] / m_Sigma;
            double besselRatio = anima::bessel_ratio_i(insideValue,m_NumberOfCoils);

            dnuValue += besselRatio * y;

            for (unsigned int j = 0;j < m_NumberOfParameters;++j)
                derivative[j] += nuDerivativeMatrix(i,j) * dnuValue;
        }
    }

}
