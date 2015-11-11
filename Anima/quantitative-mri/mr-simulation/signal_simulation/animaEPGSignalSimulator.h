#pragma once

#include <vector>

#include "AnimaSignalSimulationExport.h"

namespace anima
{
    
class ANIMASIGNALSIMULATION_EXPORT EPGSignalSimulator
{
public:
    EPGSignalSimulator();
    virtual ~EPGSignalSimulator() {}

    typedef std::vector <double> RealVectorType;

    RealVectorType GetValue(double t1Value, double t2Value,
                            double b1Value, double m0Value);

    void SetEchoSpacing(double val) {m_EchoSpacing = val;}
    void SetExcitationFlipAngle(double val) {m_ExcitationFlipAngle = val;}
    void SetFlipAngle(double val) {m_FlipAngle = val;}

    void SetNumberOfEchoes(unsigned int val) {m_NumberOfEchoes = val;}

    void SetB1OnExcitationAngle(bool val) {m_B1OnExcitationAngle = val;}

protected:
    void ComputeT2SignalMatrixElements(double t1Value, double t2Value, double b1Value);

private:
    double m_EchoSpacing;
    double m_ExcitationFlipAngle;
    double m_FlipAngle;
    unsigned int m_NumberOfEchoes;

    bool m_B1OnExcitationAngle;

    RealVectorType m_FirstLineElements, m_FirstColumnElements;
    RealVectorType m_FirstDiagonalElements, m_DiagonalElements, m_LastDiagonalElements;
    RealVectorType m_FirstLeftElements, m_FirstRightElements;
    double m_SecondLeftElement, m_SecondRightElement;
};
    
} // end namespace of anima
