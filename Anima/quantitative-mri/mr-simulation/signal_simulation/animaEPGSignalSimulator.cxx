#include <cmath>

#include <iostream>
#include "animaEPGSignalSimulator.h"

namespace anima
{
    
EPGSignalSimulator::EPGSignalSimulator()
{
    m_NumberOfEchoes = 1;
    m_EchoSpacing = 10;
    m_ExcitationFlipAngle = M_PI / 2.0;
}

EPGSignalSimulator::RealVectorType &EPGSignalSimulator::GetValue(double t1Value, double t2Value,
                                                                 double flipAngle, double m0Value)
{
    m_SimulatedT2Values.set_size(3 * m_NumberOfEchoes + 1,m_NumberOfEchoes + 1);
    m_OutputVector.resize(m_NumberOfEchoes);

    this->ComputeT2SignalMatrixElements(t1Value,t2Value,flipAngle);

    double baseValue = m0Value * std::sin(m_ExcitationFlipAngle);
    m_SimulatedT2Values(0,0) = baseValue;
    for (unsigned int i = 1;i <= m_NumberOfEchoes;++i)
        m_SimulatedT2Values(0,i) = 0.0;

    // Loop on all signals to be generated
    for (unsigned int i = 0;i < m_NumberOfEchoes;++i)
    {
        // First line
        m_SimulatedT2Values(i+1,0) = m_FirstEPGProduct * m_SimulatedT2Values(i,0) - m_SecondEPGProduct * m_SimulatedT2Values(i,3);
        if (m_NumberOfEchoes > 1)
            m_SimulatedT2Values(i+1,0) += m_ThirdEPGProduct * m_SimulatedT2Values(i,5);

        // First block
        m_SimulatedT2Values(i+1,1) = m_FourthEPGProduct * m_SimulatedT2Values(i,2);
        m_SimulatedT2Values(i+1,2) = m_FirstEPGProduct * m_SimulatedT2Values(i,1);
        m_SimulatedT2Values(i+1,3) = m_FifthEPGProduct * m_SimulatedT2Values(i,3) - m_SecondEPGProduct * m_SimulatedT2Values(i,0) / 2.0;

        if (m_NumberOfEchoes > 1)
        {
            m_SimulatedT2Values(i+1,2) -= m_SecondEPGProduct * m_SimulatedT2Values(i,6);
            m_SimulatedT2Values(i+1,3) += m_SecondEPGProduct * m_SimulatedT2Values(i,5) / 2.0;

            if (m_NumberOfEchoes > 2)
                m_SimulatedT2Values(i+1,2) += m_ThirdEPGProduct * m_SimulatedT2Values(i,8);
        }

        //center blocks
        unsigned int j;
        for (j = 1;j < m_NumberOfEchoes - 1;++j)
        {
            m_SimulatedT2Values(i+1,1 + j * 3) = m_SecondEPGProduct * m_SimulatedT2Values(i,3 + (j - 1) * 3) + m_FirstEPGProduct * m_SimulatedT2Values(i,2 + j * 3);
            if (j > 1)
                m_SimulatedT2Values(i+1,1 + j * 3) += m_ThirdEPGProduct * m_SimulatedT2Values(i,1 + (j - 2) * 3);
            else
                m_SimulatedT2Values(i+1,1 + j * 3) += m_ThirdEPGProduct * m_SimulatedT2Values(i,0);

            m_SimulatedT2Values(i+1,2 + j * 3) = m_FirstEPGProduct * m_SimulatedT2Values(i,1 + j * 3);
            m_SimulatedT2Values(i+1,3 + j * 3) = m_FifthEPGProduct * m_SimulatedT2Values(i,3 + j * 3) - m_SecondEPGProduct * m_SimulatedT2Values(i,1 + (j - 1) * 3) / 2.0;

            if ((j + 1) < m_NumberOfEchoes)
            {
                m_SimulatedT2Values(i+1,2 + j * 3) -= m_SecondEPGProduct * m_SimulatedT2Values(i,3 + (j + 1) * 3);
                m_SimulatedT2Values(i+1,3 + j * 3) += m_SecondEPGProduct * m_SimulatedT2Values(i,2 + (j + 1) * 3) / 2.0;

                if ((j + 2) < m_NumberOfEchoes)
                    m_SimulatedT2Values(i+1,2 + j * 3) += m_ThirdEPGProduct * m_SimulatedT2Values(i,2 + (j + 2) * 3);
            }
        }

        // end block line
        j = m_NumberOfEchoes - 1;
        m_SimulatedT2Values(i+1,1 + j * 3) = m_SecondEPGProduct * m_SimulatedT2Values(i,3 + (j - 1) * 3) + m_FirstEPGProduct * m_SimulatedT2Values(i,2 + j * 3);
        m_SimulatedT2Values(i+1,2 + j * 3) = 0.0;
        m_SimulatedT2Values(i+1,3 + j * 3) = m_FifthEPGProduct * m_SimulatedT2Values(i,3 + j * 3) - m_SecondEPGProduct * m_SimulatedT2Values(i,1 + (j - 1) * 3) / 2.0;

        if (j > 1)
            m_SimulatedT2Values(i+1,1 + j * 3) += m_ThirdEPGProduct * m_SimulatedT2Values(i,1 + (j - 2) * 3);
        else
            m_SimulatedT2Values(i+1,1 + j * 3) += m_ThirdEPGProduct * m_SimulatedT2Values(i,0);

        m_OutputVector[i] = m_SimulatedT2Values(i+1,0);
    }

    return m_OutputVector;
}

void EPGSignalSimulator::ComputeT2SignalMatrixElements(double t1Value, double t2Value,
                                                       double flipAngle)
{
    double espT2Value = std::exp(- m_EchoSpacing / (2 * t2Value));
    double espT1Value = std::exp(- m_EchoSpacing / (2 * t1Value));

    double cosB1alpha = std::cos(flipAngle);
    double cosB1alpha2 = std::cos(flipAngle / 2.0);
    double sinB1alpha = std::sin(flipAngle);
    double sinB1alpha2 = std::sin(flipAngle / 2.0);

    m_FirstEPGProduct = sinB1alpha2 * sinB1alpha2 * espT2Value * espT2Value;
    m_SecondEPGProduct = sinB1alpha * espT1Value * espT2Value;
    m_ThirdEPGProduct = cosB1alpha2 * cosB1alpha2 * espT2Value * espT2Value;
    m_FourthEPGProduct = espT2Value * espT2Value;
    m_FifthEPGProduct = cosB1alpha * espT1Value * espT1Value;

    m_FirstDerivativeProduct = cosB1alpha2 * sinB1alpha2 * espT2Value * espT2Value;
    m_SecondDerivativeProduct = cosB1alpha * espT1Value * espT2Value;
    m_ThirdDerivativeProduct = - sinB1alpha * espT1Value * espT1Value;
}

EPGSignalSimulator::RealVectorType &EPGSignalSimulator::GetFADerivative()
{
    m_OutputB1Derivative.resize(m_NumberOfEchoes);

    m_SimulatedDerivativeT2Values.set_size(3 * m_NumberOfEchoes + 1,m_NumberOfEchoes + 1);
    for (unsigned int i = 0;i <= m_NumberOfEchoes;++i)
        m_SimulatedDerivativeT2Values(0,i) = 0.0;

    // Loop on all signals to be generated
    for (unsigned int i = 0;i < m_NumberOfEchoes;++i)
    {
        // Start by doing dE * simulatedValues
        // First line
        m_SimulatedDerivativeT2Values(i+1,0) = m_FirstDerivativeProduct * m_SimulatedT2Values(i,0) - m_SecondDerivativeProduct * m_SimulatedT2Values(i,3);

        // First block
        m_SimulatedDerivativeT2Values(i+1,1) = 0.0;
        m_SimulatedDerivativeT2Values(i+1,2) = m_FirstDerivativeProduct * m_SimulatedT2Values(i,1);
        m_SimulatedDerivativeT2Values(i+1,3) = m_ThirdDerivativeProduct * m_SimulatedT2Values(i,3) - m_SecondDerivativeProduct * m_SimulatedT2Values(i,0) / 2.0;

        if (m_NumberOfEchoes > 1)
        {
            m_SimulatedDerivativeT2Values(i+1,0) -= m_FirstDerivativeProduct * m_SimulatedT2Values(i,5);
            m_SimulatedDerivativeT2Values(i+1,2) -= m_SecondDerivativeProduct * m_SimulatedT2Values(i,6);

            if (m_NumberOfEchoes > 2)
                m_SimulatedDerivativeT2Values(i+1,2) -= m_FirstDerivativeProduct * m_SimulatedT2Values(i,8);

            m_SimulatedDerivativeT2Values(i+1,3) += m_SecondDerivativeProduct * m_SimulatedT2Values(i,5) / 2.0;
        }

        //center blocks
        unsigned int j;
        for (j = 1;j < m_NumberOfEchoes - 1;++j)
        {
            m_SimulatedDerivativeT2Values(i+1,1 + j * 3) = m_SecondDerivativeProduct * m_SimulatedT2Values(i,3 + (j - 1) * 3) + m_FirstDerivativeProduct * m_SimulatedT2Values(i,2 + j * 3);
            if (j > 1)
                m_SimulatedDerivativeT2Values(i+1,1 + j * 3) -= m_FirstDerivativeProduct * m_SimulatedT2Values(i,1 + (j - 2) * 3);
            else
                m_SimulatedDerivativeT2Values(i+1,1 + j * 3) -= m_FirstDerivativeProduct * m_SimulatedT2Values(i,0);

            m_SimulatedDerivativeT2Values(i+1,2 + j * 3) = m_FirstDerivativeProduct * m_SimulatedT2Values(i,1 + j * 3);
            m_SimulatedDerivativeT2Values(i+1,3 + j * 3) = m_ThirdDerivativeProduct * m_SimulatedT2Values(i,3 + j * 3) - m_SecondDerivativeProduct * m_SimulatedT2Values(i,1 + (j - 1) * 3) / 2.0;

            if ((j + 1) < m_NumberOfEchoes)
            {
                m_SimulatedDerivativeT2Values(i+1,2 + j * 3) -= m_SecondDerivativeProduct * m_SimulatedT2Values(i,3 + (j + 1) * 3);
                m_SimulatedDerivativeT2Values(i+1,3 + j * 3) += m_SecondDerivativeProduct * m_SimulatedT2Values(i,2 + (j + 1) * 3) / 2.0;

                if ((j + 2) < m_NumberOfEchoes)
                    m_SimulatedDerivativeT2Values(i+1,2 + j * 3) -= m_FirstDerivativeProduct * m_SimulatedT2Values(i,2 + (j + 2) * 3);
            }
        }

        // end block line
        j = m_NumberOfEchoes - 1;
        m_SimulatedDerivativeT2Values(i+1,1 + j * 3) = m_SecondDerivativeProduct * m_SimulatedT2Values(i,3 + (j - 1) * 3) + m_FirstDerivativeProduct * m_SimulatedT2Values(i,2 + j * 3);
        if (j > 1)
            m_SimulatedDerivativeT2Values(i+1,1 + j * 3) -= m_FirstDerivativeProduct * m_SimulatedT2Values(i,1 + (j - 2) * 3);
        else
            m_SimulatedDerivativeT2Values(i+1,1 + j * 3) -= m_FirstDerivativeProduct * m_SimulatedT2Values(i,0);

        m_SimulatedDerivativeT2Values(i+1,2 + j * 3) = 0.0;
        m_SimulatedDerivativeT2Values(i+1,3 + j * 3) = m_ThirdDerivativeProduct * m_SimulatedT2Values(i,3 + j * 3) - m_SecondDerivativeProduct * m_SimulatedT2Values(i,1 + (j - 1) * 3) / 2.0;

        // Now adding E * simulatedDerivative[i]
        // First line
        m_SimulatedDerivativeT2Values(i+1,0) += m_FirstEPGProduct * m_SimulatedDerivativeT2Values(i,0) - m_SecondEPGProduct * m_SimulatedDerivativeT2Values(i,3);
        if (m_NumberOfEchoes > 1)
            m_SimulatedDerivativeT2Values(i+1,0) += m_ThirdEPGProduct * m_SimulatedDerivativeT2Values(i,5);

        // First block
        m_SimulatedDerivativeT2Values(i+1,1) += m_FourthEPGProduct * m_SimulatedDerivativeT2Values(i,2);
        m_SimulatedDerivativeT2Values(i+1,2) += m_FirstEPGProduct * m_SimulatedDerivativeT2Values(i,1);
        m_SimulatedDerivativeT2Values(i+1,3) += m_FifthEPGProduct * m_SimulatedDerivativeT2Values(i,3) - m_SecondEPGProduct * m_SimulatedDerivativeT2Values(i,0) / 2.0;
        if (m_NumberOfEchoes > 1)
        {
            m_SimulatedDerivativeT2Values(i+1,2) -= m_SecondEPGProduct * m_SimulatedDerivativeT2Values(i,6);
            m_SimulatedDerivativeT2Values(i+1,3) += m_SecondEPGProduct * m_SimulatedDerivativeT2Values(i,5) / 2.0;

            if (m_NumberOfEchoes > 2)
                m_SimulatedDerivativeT2Values(i+1,2) += m_ThirdEPGProduct * m_SimulatedDerivativeT2Values(i,8);
        }

        //center blocks
        for (j = 1;j < m_NumberOfEchoes - 1;++j)
        {
            m_SimulatedDerivativeT2Values(i+1,1 + j * 3) += m_SecondEPGProduct * m_SimulatedDerivativeT2Values(i,3 + (j - 1) * 3) + m_FirstEPGProduct * m_SimulatedDerivativeT2Values(i,2 + j * 3);
            if (j > 1)
                m_SimulatedDerivativeT2Values(i+1,1 + j * 3) += m_ThirdEPGProduct * m_SimulatedDerivativeT2Values(i,1 + (j - 2) * 3);
            else
                m_SimulatedDerivativeT2Values(i+1,1 + j * 3) += m_ThirdEPGProduct * m_SimulatedDerivativeT2Values(i,0);

            m_SimulatedDerivativeT2Values(i+1,2 + j * 3) += m_FirstEPGProduct * m_SimulatedDerivativeT2Values(i,1 + j * 3);
            m_SimulatedDerivativeT2Values(i+1,3 + j * 3) += m_FifthEPGProduct * m_SimulatedDerivativeT2Values(i,3 + j * 3) - m_SecondEPGProduct * m_SimulatedDerivativeT2Values(i,1 + (j - 1) * 3) / 2.0;

            if ((j + 1) < m_NumberOfEchoes)
            {
                m_SimulatedDerivativeT2Values(i+1,2 + j * 3) -= m_SecondEPGProduct * m_SimulatedDerivativeT2Values(i,3 + (j + 1) * 3);
                m_SimulatedDerivativeT2Values(i+1,3 + j * 3) += m_SecondEPGProduct * m_SimulatedDerivativeT2Values(i,2 + (j + 1) * 3) / 2.0;

                if ((j + 2) < m_NumberOfEchoes)
                    m_SimulatedDerivativeT2Values(i+1,2 + j * 3) += m_ThirdEPGProduct * m_SimulatedDerivativeT2Values(i,2 + (j + 2) * 3);
            }
        }

        // end block line
        j = m_NumberOfEchoes - 1;
        m_SimulatedDerivativeT2Values(i+1,1 + j * 3) += m_SecondEPGProduct * m_SimulatedDerivativeT2Values(i,3 + (j - 1) * 3) + m_FirstEPGProduct * m_SimulatedDerivativeT2Values(i,2 + j * 3);
        m_SimulatedDerivativeT2Values(i+1,3 + j * 3) += m_FifthEPGProduct * m_SimulatedDerivativeT2Values(i,3 + j * 3) - m_SecondEPGProduct * m_SimulatedDerivativeT2Values(i,1 + (j - 1) * 3) / 2.0;

        if (j > 1)
            m_SimulatedDerivativeT2Values(i+1,1 + j * 3) += m_ThirdEPGProduct * m_SimulatedDerivativeT2Values(i,1 + (j - 2) * 3);
        else
            m_SimulatedDerivativeT2Values(i+1,1 + j * 3) += m_ThirdEPGProduct * m_SimulatedDerivativeT2Values(i,0);

        m_OutputB1Derivative[i] = m_SimulatedDerivativeT2Values(i+1,0);
    }

    return m_OutputB1Derivative;
}
    
} // end of namespace anima
