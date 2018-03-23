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
    m_FlipAngle = M_PI;
    m_B1OnExcitationAngle = false;

    m_FirstLineElements.resize(3);
    m_FirstColumnElements.resize(2);
    m_FirstDiagonalElements.resize(3);
    m_DiagonalElements.resize(3);
    m_LastDiagonalElements.resize(2);
    m_FirstLeftElements.resize(2);
    m_FirstRightElements.resize(2);
}

EPGSignalSimulator::RealVectorType &EPGSignalSimulator::GetValue(double t1Value, double t2Value,
                                                                double b1Value, double m0Value)
{
    m_SimulatedT2Values.resize(3 * m_NumberOfEchoes + 1);
    m_WorkVector.resize(3 * m_NumberOfEchoes + 1);
    m_OutputVector.resize(m_NumberOfEchoes);

    std::fill(m_SimulatedT2Values.begin(),m_SimulatedT2Values.end(),0.0);
    std::fill(m_OutputVector.begin(),m_OutputVector.end(),0.0);
    std::fill(m_WorkVector.begin(),m_WorkVector.end(),0.0);

    this->ComputeT2SignalMatrixElements(t1Value,t2Value,b1Value);

    double baseValue;
    if (m_B1OnExcitationAngle)
        baseValue = m0Value * std::sin(b1Value * m_ExcitationFlipAngle);
    else
        baseValue = m0Value * std::sin(m_ExcitationFlipAngle);

    m_WorkVector[0] = baseValue;

    // Loop on all signals to be generated
    for (unsigned int i = 0;i < m_NumberOfEchoes;++i)
    {
        std::copy(m_WorkVector.begin(),m_WorkVector.end(),m_SimulatedT2Values.begin());
        std::fill(m_WorkVector.begin(),m_WorkVector.end(),0.0);

        // First line
        m_WorkVector[0] = m_FirstLineElements[0] * m_SimulatedT2Values[0] + m_FirstLineElements[1] * m_SimulatedT2Values[3];
        if (m_NumberOfEchoes > 1)
            m_WorkVector[0] += m_FirstLineElements[2] * m_SimulatedT2Values[5];

        // First block
        m_WorkVector[1] = m_FirstDiagonalElements[0] * m_SimulatedT2Values[2];
        m_WorkVector[2] = m_FirstDiagonalElements[1] * m_SimulatedT2Values[1];
        if (m_NumberOfEchoes > 1)
            m_WorkVector[2] += m_FirstRightElements[0] * m_SimulatedT2Values[6];
        if (m_NumberOfEchoes > 2)
            m_WorkVector[2] += m_SecondRightElement * m_SimulatedT2Values[8];
        m_WorkVector[3] = m_FirstColumnElements[0] * m_SimulatedT2Values[0] + m_FirstDiagonalElements[2] * m_SimulatedT2Values[3];
        if (m_NumberOfEchoes > 1)
            m_WorkVector[3] += m_FirstRightElements[1] * m_SimulatedT2Values[5];

        //center blocks
        for (unsigned int j = 1;j < m_NumberOfEchoes - 1;++j)
        {
            m_WorkVector[1 + j * 3] = m_FirstLeftElements[0] * m_SimulatedT2Values[3 + (j - 1) * 3] + m_DiagonalElements[0] * m_SimulatedT2Values[2 + j * 3];
            if (j > 1)
                m_WorkVector[1 + j * 3] += m_SecondLeftElement * m_SimulatedT2Values[1 + (j - 2) * 3];
            else
                m_WorkVector[1 + j * 3] += m_FirstColumnElements[1] * m_SimulatedT2Values[0];

            m_WorkVector[2 + j * 3] = m_DiagonalElements[1] * m_SimulatedT2Values[1 + j * 3];

            if ((j + 1) < m_NumberOfEchoes)
                m_WorkVector [2 + j * 3] += m_FirstRightElements[0] * m_SimulatedT2Values[3 + (j + 1) * 3];

            if ((j + 2) < m_NumberOfEchoes)
                m_WorkVector [2 + j * 3] += m_SecondRightElement * m_SimulatedT2Values[2 + (j + 2) * 3];

            m_WorkVector[3 + j * 3] = m_FirstLeftElements[1] * m_SimulatedT2Values[1 + (j - 1) * 3] + m_DiagonalElements[2] * m_SimulatedT2Values[3 + j * 3];

            if ((j + 1) < m_NumberOfEchoes)
                m_WorkVector [3 + j * 3] += m_FirstRightElements[1] * m_SimulatedT2Values[2 + (j + 1) * 3];
        }

        // end block line
        unsigned int j = m_NumberOfEchoes - 1;
        m_WorkVector[1 + j * 3] = m_FirstLeftElements[0] * m_SimulatedT2Values[3 + (j - 1) * 3] + m_LastDiagonalElements[0] * m_SimulatedT2Values[2 + j * 3];
        if (j > 1)
            m_WorkVector[1 + j * 3] += m_SecondLeftElement * m_SimulatedT2Values[1 + (j - 2) * 3];
        else
            m_WorkVector[1 + j * 3] += m_FirstColumnElements[1] * m_SimulatedT2Values[0];

        m_WorkVector[2 + j * 3] = 0;
        m_WorkVector[3 + j * 3] = m_FirstLeftElements[1] * m_SimulatedT2Values[1 + (j - 1) * 3] + m_LastDiagonalElements[1] * m_SimulatedT2Values[3 + j * 3];

        m_OutputVector[i] = m_WorkVector[0];
    }

    return m_OutputVector;
}

void EPGSignalSimulator::ComputeT2SignalMatrixElements(double t1Value, double t2Value,
                                                       double b1Value)
{
    double espT2Value = std::exp(- m_EchoSpacing / (2 * t2Value));
    double espT1Value = std::exp(- m_EchoSpacing / (2 * t1Value));

    double cosB1alpha = std::cos(b1Value * m_FlipAngle);
    double cosB1alpha2 = std::cos(b1Value * m_FlipAngle / 2.0);
    double sinB1alpha = std::sin(b1Value * m_FlipAngle);
    double sinB1alpha2 = std::sin(b1Value * m_FlipAngle / 2.0);

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
    
} // end of namespace anima
