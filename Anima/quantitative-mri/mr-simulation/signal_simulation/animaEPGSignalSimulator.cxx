#define _USE_MATH_DEFINES
#include <cmath>

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
}

EPGSignalSimulator::RealVectorType EPGSignalSimulator::GetValue(double t1Value, double t2Value,
                                                                double b1Value, double m0Value)
{
    RealVectorType outputVector(m_NumberOfEchoes,0);

    this->ComputeT2SignalMatrices(t1Value,t2Value,b1Value);

    std::complex <double> zeroValue(0,0);
    ComplexVectorType simulatedT2Values(3 * m_NumberOfEchoes,zeroValue);
    ComplexVectorType workT2Values(3 * m_NumberOfEchoes,zeroValue);

    double baseValue;
    if (m_B1OnExcitationAngle)
        baseValue = m0Value * std::sin(b1Value * m_ExcitationFlipAngle);
    else
        baseValue = m0Value * std::sin(m_ExcitationFlipAngle);

    simulatedT2Values[0] = std::complex<double>(baseValue,0);

    // Loop on all signals to be generated
    for (unsigned int i = 0;i < m_NumberOfEchoes;++i)
    {
        unsigned int maxJIndex = i + 1;
        if (maxJIndex > m_NumberOfEchoes - 1)
            maxJIndex = m_NumberOfEchoes - 1;

        for (unsigned int j = 0;j <= maxJIndex;++j)
        {
            // Multiply j-1 sub-matrix by lower diagonal term
            if (j > 0)
            {
                for (unsigned int k = 0;k < 3;++k)
                    workT2Values[3 * j] += m_LowerDiagonalT2SignalMatrices[i][k] * simulatedT2Values[3 * (j-1) + k];
            }

            // Multiply j sub-matrix by diagonal term
            if (j == 0)
            {
                for (unsigned int k = 0;k < 3;++k)
                    workT2Values[3*j] += m_DiagonalT2SignalMatrices[i][0][k] * simulatedT2Values[3*j + k];
            }

            for (unsigned int k = 0;k < 3;++k)
                workT2Values[3*j + 2] += m_DiagonalT2SignalMatrices[i][2][k] * simulatedT2Values[3*j + k];

            // Multiply j+1 sub-matrix by upper diagonal term
            if (j < (m_NumberOfEchoes - 1))
            {
                for (unsigned int k = 0;k < 3;++k)
                    workT2Values[3*j + 1] += m_UpperDiagonalT2SignalMatrices[i][k] * simulatedT2Values[3*(j+1) + k];
            }
        }

        simulatedT2Values = workT2Values;
        std::fill(workT2Values.begin(),workT2Values.end(),zeroValue);

        outputVector[i] = simulatedT2Values[0].real();
    }

    return outputVector;
}

void EPGSignalSimulator::ComputeT2SignalMatrices(double t1Value, double t2Value,
                                                 double b1Value)
{
    m_DiagonalT2SignalMatrices.resize(m_NumberOfEchoes);
    m_LowerDiagonalT2SignalMatrices.resize(m_NumberOfEchoes);
    m_UpperDiagonalT2SignalMatrices.resize(m_NumberOfEchoes);

    unsigned int matrixSize = 3;

    double espT2Value = std::exp(- m_EchoSpacing / (2 * t2Value));
    double espT1Value = std::exp(- m_EchoSpacing / (2 * t1Value));

    std::complex <double> zeroValue(0,0);
    ComplexVectorType zeroVector(matrixSize,zeroValue);

    for (unsigned int i = 0;i < m_NumberOfEchoes;++i)
    {
        m_DiagonalT2SignalMatrices[i].resize(matrixSize);

        for (unsigned int j = 0;j < matrixSize;++j)
            m_DiagonalT2SignalMatrices[i][j] = zeroVector;

        double cosB1alpha = std::cos(b1Value * m_FlipAngle);
        double cosB1alpha2 = std::cos(b1Value * m_FlipAngle / 2.0);
        double sinB1alpha = std::sin(b1Value * m_FlipAngle);
        double sinB1alpha2 = std::sin(b1Value * m_FlipAngle / 2.0);

        m_DiagonalT2SignalMatrices[i][0][0] = std::complex<double>(espT2Value * espT2Value * sinB1alpha2 * sinB1alpha2,0);
        m_DiagonalT2SignalMatrices[i][0][1] = std::complex<double>(espT2Value * espT2Value * cosB1alpha2 * cosB1alpha2,0);
        m_DiagonalT2SignalMatrices[i][0][2] = std::complex<double>(0,espT2Value * espT1Value * sinB1alpha);

        m_DiagonalT2SignalMatrices[i][2][0] = std::complex<double>(0,- 0.5 * espT2Value * espT1Value * sinB1alpha);
        m_DiagonalT2SignalMatrices[i][2][1] = std::complex<double>(0,0.5 * espT2Value * espT1Value * sinB1alpha);
        m_DiagonalT2SignalMatrices[i][2][2] = std::complex<double>(espT1Value * espT1Value * cosB1alpha,0);

        m_LowerDiagonalT2SignalMatrices[i] = zeroVector;
        m_LowerDiagonalT2SignalMatrices[i][0] = std::complex<double>(espT2Value * espT2Value * cosB1alpha2 * cosB1alpha2,0);
        m_LowerDiagonalT2SignalMatrices[i][1] = std::complex<double>(espT2Value * espT2Value * sinB1alpha2 * sinB1alpha2,0);
        m_LowerDiagonalT2SignalMatrices[i][2] = std::complex<double>(0,- espT2Value * espT1Value * sinB1alpha);

        m_UpperDiagonalT2SignalMatrices[i] = zeroVector;
        m_UpperDiagonalT2SignalMatrices[i][0] = std::complex<double>(espT2Value * espT2Value * sinB1alpha2 * sinB1alpha2,0);
        m_UpperDiagonalT2SignalMatrices[i][1] = std::complex<double>(espT2Value * espT2Value * cosB1alpha2 * cosB1alpha2,0);
        m_UpperDiagonalT2SignalMatrices[i][2] = std::complex<double>(0,espT2Value * espT1Value * sinB1alpha);
    }
}
    
} // end of namespace anima
