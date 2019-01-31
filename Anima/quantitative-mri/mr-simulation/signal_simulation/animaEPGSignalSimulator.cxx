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

    m_FirstLineElements.resize(3);
    m_FirstColumnElements.resize(2);
    m_FirstDiagonalElements.resize(3);
    m_DiagonalElements.resize(3);
    m_LastDiagonalElements.resize(2);
    m_FirstLeftElements.resize(2);
    m_FirstRightElements.resize(2);

    m_FirstLineDerivativeElements.resize(3);
    m_FirstColumnDerivativeElements.resize(2);
    m_FirstDiagonalDerivativeElements.resize(2);
    m_DiagonalDerivativeElements.resize(3);
    m_LastDiagonalDerivativeElements.resize(2);
    m_FirstLeftDerivativeElements.resize(2);
    m_FirstRightDerivativeElements.resize(2);
}

EPGSignalSimulator::RealVectorType &EPGSignalSimulator::GetValue(double t1Value, double t2Value,
                                                                double flipAngle, double m0Value)
{
    m_SimulatedT2Values.set_size(3 * m_NumberOfEchoes + 1,m_NumberOfEchoes + 1);
    m_OutputVector.resize(m_NumberOfEchoes);

    m_SimulatedT2Values.fill(0.0);
    std::fill(m_OutputVector.begin(),m_OutputVector.end(),0.0);

    this->ComputeT2SignalMatrixElements(t1Value,t2Value,flipAngle);

    double baseValue = m0Value * std::sin(m_ExcitationFlipAngle);
    m_SimulatedT2Values(0,0) = baseValue;

    // Loop on all signals to be generated
    for (unsigned int i = 0;i < m_NumberOfEchoes;++i)
    {
        // First line
        m_SimulatedT2Values(i+1,0) = m_FirstLineElements[0] * m_SimulatedT2Values(i,0) + m_FirstLineElements[1] * m_SimulatedT2Values(i,3);
        if (m_NumberOfEchoes > 1)
            m_SimulatedT2Values(i+1,0) += m_FirstLineElements[2] * m_SimulatedT2Values(i,5);

        // First block
        m_SimulatedT2Values(i+1,1) = m_FirstDiagonalElements[0] * m_SimulatedT2Values(i,2);
        m_SimulatedT2Values(i+1,2) = m_FirstDiagonalElements[1] * m_SimulatedT2Values(i,1);
        if (m_NumberOfEchoes > 1)
            m_SimulatedT2Values(i+1,2) += m_FirstRightElements[0] * m_SimulatedT2Values(i,6);
        if (m_NumberOfEchoes > 2)
            m_SimulatedT2Values(i+1,2) += m_SecondRightElement * m_SimulatedT2Values(i,8);
        m_SimulatedT2Values(i+1,3) = m_FirstColumnElements[0] * m_SimulatedT2Values(i,0) + m_FirstDiagonalElements[2] * m_SimulatedT2Values(i,3);
        if (m_NumberOfEchoes > 1)
            m_SimulatedT2Values(i+1,3) += m_FirstRightElements[1] * m_SimulatedT2Values(i,5);

        //center blocks
        for (unsigned int j = 1;j < m_NumberOfEchoes - 1;++j)
        {
            m_SimulatedT2Values(i+1,1 + j * 3) = m_FirstLeftElements[0] * m_SimulatedT2Values(i,3 + (j - 1) * 3) + m_DiagonalElements[0] * m_SimulatedT2Values(i,2 + j * 3);
            if (j > 1)
                m_SimulatedT2Values(i+1,1 + j * 3) += m_SecondLeftElement * m_SimulatedT2Values(i,1 + (j - 2) * 3);
            else
                m_SimulatedT2Values(i+1,1 + j * 3) += m_FirstColumnElements[1] * m_SimulatedT2Values(i,0);

            m_SimulatedT2Values(i+1,2 + j * 3) = m_DiagonalElements[1] * m_SimulatedT2Values(i,1 + j * 3);

            if ((j + 1) < m_NumberOfEchoes)
                m_SimulatedT2Values(i+1,2 + j * 3) += m_FirstRightElements[0] * m_SimulatedT2Values(i,3 + (j + 1) * 3);

            if ((j + 2) < m_NumberOfEchoes)
                m_SimulatedT2Values(i+1,2 + j * 3) += m_SecondRightElement * m_SimulatedT2Values(i,2 + (j + 2) * 3);

            m_SimulatedT2Values(i+1,3 + j * 3) = m_FirstLeftElements[1] * m_SimulatedT2Values(i,1 + (j - 1) * 3) + m_DiagonalElements[2] * m_SimulatedT2Values(i,3 + j * 3);

            if ((j + 1) < m_NumberOfEchoes)
                m_SimulatedT2Values(i+1,3 + j * 3) += m_FirstRightElements[1] * m_SimulatedT2Values(i,2 + (j + 1) * 3);
        }

        // end block line
        unsigned int j = m_NumberOfEchoes - 1;
        m_SimulatedT2Values(i+1,1 + j * 3) = m_FirstLeftElements[0] * m_SimulatedT2Values(i,3 + (j - 1) * 3) + m_LastDiagonalElements[0] * m_SimulatedT2Values(i,2 + j * 3);
        if (j > 1)
            m_SimulatedT2Values(i+1,1 + j * 3) += m_SecondLeftElement * m_SimulatedT2Values(i,1 + (j - 2) * 3);
        else
            m_SimulatedT2Values(i+1,1 + j * 3) += m_FirstColumnElements[1] * m_SimulatedT2Values(i,0);

        m_SimulatedT2Values(i+1,3 + j * 3) = m_FirstLeftElements[1] * m_SimulatedT2Values(i,1 + (j - 1) * 3) + m_LastDiagonalElements[1] * m_SimulatedT2Values(i,3 + j * 3);

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

    m_FirstLineElements[0] = sinB1alpha2 * sinB1alpha2 * espT2Value * espT2Value;
    m_FirstLineElements[1] = - sinB1alpha * espT1Value * espT2Value;
    m_FirstLineElements[2] = cosB1alpha2 * cosB1alpha2 * espT2Value * espT2Value;

    m_FirstColumnElements[0] = - sinB1alpha * espT1Value * espT2Value / 2.0;
    m_FirstColumnElements[1] = cosB1alpha2 * cosB1alpha2 * espT2Value * espT2Value;

    m_FirstDiagonalElements[0] = espT2Value * espT2Value;
    m_FirstDiagonalElements[1] = sinB1alpha2 * sinB1alpha2 * espT2Value * espT2Value;
    m_FirstDiagonalElements[2] = cosB1alpha * espT1Value * espT1Value;

    m_DiagonalElements[0] = sinB1alpha2 * sinB1alpha2 * espT2Value * espT2Value;
    m_DiagonalElements[1] = sinB1alpha2 * sinB1alpha2 * espT2Value * espT2Value;
    m_DiagonalElements[2] = cosB1alpha * espT1Value * espT1Value;

    m_LastDiagonalElements[0] = sinB1alpha2 * sinB1alpha2 * espT2Value * espT2Value;
    m_LastDiagonalElements[1] = cosB1alpha * espT1Value * espT1Value;

    m_FirstRightElements[0] = - sinB1alpha * espT1Value * espT2Value;
    m_FirstRightElements[1] = sinB1alpha * espT1Value * espT2Value / 2.0;

    m_SecondRightElement = cosB1alpha2 * cosB1alpha2 * espT2Value * espT2Value;

    m_FirstLeftElements[0] = sinB1alpha * espT1Value * espT2Value;
    m_FirstLeftElements[1] = - sinB1alpha * espT1Value * espT2Value / 2.0;

    m_SecondLeftElement = cosB1alpha2 * cosB1alpha2 * espT2Value * espT2Value;
}

EPGSignalSimulator::RealVectorType &EPGSignalSimulator::GetFADerivative(double t1Value, double t2Value,
                                                                        double flipAngle, double m0Value)
{
    m_OutputB1Derivative.resize(m_NumberOfEchoes);
    std::fill(m_OutputB1Derivative.begin(),m_OutputB1Derivative.end(),0.0);

    m_SimulatedDerivativeT2Values.set_size(3 * m_NumberOfEchoes + 1,m_NumberOfEchoes + 1);
    m_SimulatedDerivativeT2Values.fill(0.0);

    this->ComputeT2SignalFADerivativeMatrixElements(t1Value,t2Value,flipAngle);

    // Loop on all signals to be generated
    for (unsigned int i = 0;i < m_NumberOfEchoes;++i)
    {
        // Start by doing dE * simulatedValues

        // First line
        m_SimulatedDerivativeT2Values(i+1,0) = m_FirstLineDerivativeElements[0] * m_SimulatedT2Values(i,0) + m_FirstLineDerivativeElements[1] * m_SimulatedT2Values(i,3);
        if (m_NumberOfEchoes > 1)
            m_SimulatedDerivativeT2Values(i+1,0) += m_FirstLineDerivativeElements[2] * m_SimulatedT2Values(i,5);

        // First block
        m_SimulatedDerivativeT2Values(i+1,2) = m_FirstDiagonalDerivativeElements[0] * m_SimulatedT2Values(i,1);
        if (m_NumberOfEchoes > 1)
            m_SimulatedDerivativeT2Values(i+1,2) += m_FirstRightDerivativeElements[0] * m_SimulatedT2Values(i,6);
        if (m_NumberOfEchoes > 2)
            m_SimulatedDerivativeT2Values(i+1,2) += m_SecondRightDerivativeElement * m_SimulatedT2Values(i,8);
        m_SimulatedDerivativeT2Values(i+1,3) = m_FirstColumnDerivativeElements[0] * m_SimulatedT2Values(i,0) + m_FirstDiagonalDerivativeElements[1] * m_SimulatedT2Values(i,3);
        if (m_NumberOfEchoes > 1)
            m_SimulatedDerivativeT2Values(i+1,3) += m_FirstRightDerivativeElements[1] * m_SimulatedT2Values(i,5);

        //center blocks
        for (unsigned int j = 1;j < m_NumberOfEchoes - 1;++j)
        {
            m_SimulatedDerivativeT2Values(i+1,1 + j * 3) = m_FirstLeftDerivativeElements[0] * m_SimulatedT2Values(i,3 + (j - 1) * 3) + m_DiagonalDerivativeElements[0] * m_SimulatedT2Values(i,2 + j * 3);
            if (j > 1)
                m_SimulatedDerivativeT2Values(i+1,1 + j * 3) += m_SecondLeftDerivativeElement * m_SimulatedT2Values(i,1 + (j - 2) * 3);
            else
                m_SimulatedDerivativeT2Values(i+1,1 + j * 3) += m_FirstColumnDerivativeElements[1] * m_SimulatedT2Values(i,0);

            m_SimulatedDerivativeT2Values(i+1,2 + j * 3) = m_DiagonalDerivativeElements[1] * m_SimulatedT2Values(i,1 + j * 3);

            if ((j + 1) < m_NumberOfEchoes)
                m_SimulatedDerivativeT2Values(i+1,2 + j * 3) += m_FirstRightDerivativeElements[0] * m_SimulatedT2Values(i,3 + (j + 1) * 3);

            if ((j + 2) < m_NumberOfEchoes)
                m_SimulatedDerivativeT2Values(i+1,2 + j * 3) += m_SecondRightDerivativeElement * m_SimulatedT2Values(i,2 + (j + 2) * 3);

            m_SimulatedDerivativeT2Values(i+1,3 + j * 3) = m_FirstLeftDerivativeElements[1] * m_SimulatedT2Values(i,1 + (j - 1) * 3) + m_DiagonalDerivativeElements[2] * m_SimulatedT2Values(i,3 + j * 3);

            if ((j + 1) < m_NumberOfEchoes)
                m_SimulatedDerivativeT2Values(i+1,3 + j * 3) += m_FirstRightDerivativeElements[1] * m_SimulatedT2Values(i,2 + (j + 1) * 3);
        }

        // end block line
        unsigned int j = m_NumberOfEchoes - 1;
        m_SimulatedDerivativeT2Values(i+1,1 + j * 3) = m_FirstLeftDerivativeElements[0] * m_SimulatedT2Values(i,3 + (j - 1) * 3) + m_LastDiagonalDerivativeElements[0] * m_SimulatedT2Values(i,2 + j * 3);
        if (j > 1)
            m_SimulatedDerivativeT2Values(i+1,1 + j * 3) += m_SecondLeftDerivativeElement * m_SimulatedT2Values(i,1 + (j - 2) * 3);
        else
            m_SimulatedDerivativeT2Values(i+1,1 + j * 3) += m_FirstColumnDerivativeElements[1] * m_SimulatedT2Values(i,0);

        m_SimulatedDerivativeT2Values(i+1,3 + j * 3) = m_FirstLeftDerivativeElements[1] * m_SimulatedT2Values(i,1 + (j - 1) * 3) + m_LastDiagonalDerivativeElements[1] * m_SimulatedT2Values(i,3 + j * 3);

        // Now adding E * simulatedDerivative[i]

        // First line
        m_SimulatedDerivativeT2Values(i+1,0) += m_FirstLineElements[0] * m_SimulatedDerivativeT2Values(i,0) + m_FirstLineElements[1] * m_SimulatedDerivativeT2Values(i,3);
        if (m_NumberOfEchoes > 1)
            m_SimulatedDerivativeT2Values(i+1,0) += m_FirstLineElements[2] * m_SimulatedDerivativeT2Values(i,5);

        // First block
        m_SimulatedDerivativeT2Values(i+1,1) += m_FirstDiagonalElements[0] * m_SimulatedDerivativeT2Values(i,2);
        m_SimulatedDerivativeT2Values(i+1,2) += m_FirstDiagonalElements[1] * m_SimulatedDerivativeT2Values(i,1);
        if (m_NumberOfEchoes > 1)
            m_SimulatedDerivativeT2Values(i+1,2) += m_FirstRightElements[0] * m_SimulatedDerivativeT2Values(i,6);
        if (m_NumberOfEchoes > 2)
            m_SimulatedDerivativeT2Values(i+1,2) += m_SecondRightElement * m_SimulatedDerivativeT2Values(i,8);
        m_SimulatedDerivativeT2Values(i+1,3) += m_FirstColumnElements[0] * m_SimulatedDerivativeT2Values(i,0) + m_FirstDiagonalElements[2] * m_SimulatedDerivativeT2Values(i,3);
        if (m_NumberOfEchoes > 1)
            m_SimulatedDerivativeT2Values(i+1,3) += m_FirstRightElements[1] * m_SimulatedDerivativeT2Values(i,5);

        //center blocks
        for (unsigned int j = 1;j < m_NumberOfEchoes - 1;++j)
        {
            m_SimulatedDerivativeT2Values(i+1,1 + j * 3) += m_FirstLeftElements[0] * m_SimulatedDerivativeT2Values(i,3 + (j - 1) * 3) + m_DiagonalElements[0] * m_SimulatedDerivativeT2Values(i,2 + j * 3);
            if (j > 1)
                m_SimulatedDerivativeT2Values(i+1,1 + j * 3) += m_SecondLeftElement * m_SimulatedDerivativeT2Values(i,1 + (j - 2) * 3);
            else
                m_SimulatedDerivativeT2Values(i+1,1 + j * 3) += m_FirstColumnElements[1] * m_SimulatedDerivativeT2Values(i,0);

            m_SimulatedDerivativeT2Values(i+1,2 + j * 3) += m_DiagonalElements[1] * m_SimulatedDerivativeT2Values(i,1 + j * 3);

            if ((j + 1) < m_NumberOfEchoes)
                m_SimulatedDerivativeT2Values(i+1,2 + j * 3) += m_FirstRightElements[0] * m_SimulatedDerivativeT2Values(i,3 + (j + 1) * 3);

            if ((j + 2) < m_NumberOfEchoes)
                m_SimulatedDerivativeT2Values(i+1,2 + j * 3) += m_SecondRightElement * m_SimulatedDerivativeT2Values(i,2 + (j + 2) * 3);

            m_SimulatedDerivativeT2Values(i+1,3 + j * 3) += m_FirstLeftElements[1] * m_SimulatedDerivativeT2Values(i,1 + (j - 1) * 3) + m_DiagonalElements[2] * m_SimulatedDerivativeT2Values(i,3 + j * 3);

            if ((j + 1) < m_NumberOfEchoes)
                m_SimulatedDerivativeT2Values(i+1,3 + j * 3) += m_FirstRightElements[1] * m_SimulatedDerivativeT2Values(i,2 + (j + 1) * 3);
        }

        // end block line
        j = m_NumberOfEchoes - 1;
        m_SimulatedDerivativeT2Values(i+1,1 + j * 3) += m_FirstLeftElements[0] * m_SimulatedDerivativeT2Values(i,3 + (j - 1) * 3) + m_LastDiagonalElements[0] * m_SimulatedDerivativeT2Values(i,2 + j * 3);
        if (j > 1)
            m_SimulatedDerivativeT2Values(i+1,1 + j * 3) += m_SecondLeftElement * m_SimulatedDerivativeT2Values(i,1 + (j - 2) * 3);
        else
            m_SimulatedDerivativeT2Values(i+1,1 + j * 3) += m_FirstColumnElements[1] * m_SimulatedDerivativeT2Values(i,0);

        m_SimulatedDerivativeT2Values(i+1,3 + j * 3) += m_FirstLeftElements[1] * m_SimulatedDerivativeT2Values(i,1 + (j - 1) * 3) + m_LastDiagonalElements[1] * m_SimulatedDerivativeT2Values(i,3 + j * 3);

        m_OutputB1Derivative[i] = m_SimulatedDerivativeT2Values(i+1,0);
    }

    return m_OutputB1Derivative;
}

void EPGSignalSimulator::ComputeT2SignalFADerivativeMatrixElements(double t1Value, double t2Value,
                                                                   double flipAngle)
{
    double espT2Value = std::exp(- m_EchoSpacing / (2 * t2Value));
    double espT1Value = std::exp(- m_EchoSpacing / (2 * t1Value));

    double cosB1alpha = std::cos(flipAngle);
    double cosB1alpha2 = std::cos(flipAngle / 2.0);
    double sinB1alpha = std::sin(flipAngle);
    double sinB1alpha2 = std::sin(flipAngle / 2.0);

    m_FirstLineDerivativeElements[0] = cosB1alpha2 * sinB1alpha2 * espT2Value * espT2Value;
    m_FirstLineDerivativeElements[1] = - cosB1alpha * espT1Value * espT2Value;
    m_FirstLineDerivativeElements[2] = - m_FirstLineDerivativeElements[0];

    m_FirstColumnDerivativeElements[0] = m_FirstLineDerivativeElements[1] / 2.0;
    m_FirstColumnDerivativeElements[1] = m_FirstLineDerivativeElements[2];

    m_FirstDiagonalDerivativeElements[0] = m_FirstLineDerivativeElements[0];
    m_FirstDiagonalDerivativeElements[1] = - sinB1alpha * espT1Value * espT1Value;

    m_DiagonalDerivativeElements[0] = m_FirstLineDerivativeElements[0];
    m_DiagonalDerivativeElements[1] = m_FirstLineDerivativeElements[0];
    m_DiagonalDerivativeElements[2] = m_FirstDiagonalDerivativeElements[1];

    m_LastDiagonalDerivativeElements[0] = m_FirstLineDerivativeElements[0];
    m_LastDiagonalDerivativeElements[1] = m_FirstDiagonalDerivativeElements[1];

    m_FirstRightDerivativeElements[0] = m_FirstLineDerivativeElements[1];
    m_FirstRightDerivativeElements[1] = - m_FirstLineDerivativeElements[1] / 2.0;

    m_SecondRightDerivativeElement = - m_FirstLineDerivativeElements[0];

    m_FirstLeftDerivativeElements[0] = - m_FirstLineDerivativeElements[1];
    m_FirstLeftDerivativeElements[1] = m_FirstLineDerivativeElements[1] / 2.0;

    m_SecondLeftDerivativeElement = - m_FirstLineDerivativeElements[0];
}
    
} // end of namespace anima
