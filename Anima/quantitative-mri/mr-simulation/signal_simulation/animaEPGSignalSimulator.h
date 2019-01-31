#pragma once

#include <vector>
#include <vnl/vnl_matrix.h>

#include "AnimaSignalSimulationExport.h"

namespace anima
{
    
class ANIMASIGNALSIMULATION_EXPORT EPGSignalSimulator
{
public:
    EPGSignalSimulator();
    virtual ~EPGSignalSimulator() {}

    typedef std::vector <double> RealVectorType;

    RealVectorType &GetValue(double t1Value, double t2Value,
                             double flipAngle, double m0Value);

    RealVectorType &GetFADerivative(double t1Value, double t2Value,
                                    double flipAngle, double m0Value);

    void SetEchoSpacing(double val) {m_EchoSpacing = val;}
    void SetExcitationFlipAngle(double val) {m_ExcitationFlipAngle = val;}

    void SetNumberOfEchoes(unsigned int val) {m_NumberOfEchoes = val;}

protected:
    void ComputeT2SignalMatrixElements(double t1Value, double t2Value, double flipAngle);
    void ComputeT2SignalFADerivativeMatrixElements(double t1Value, double t2Value, double flipAngle);

private:
    double m_EchoSpacing;
    double m_ExcitationFlipAngle;
    unsigned int m_NumberOfEchoes;

    RealVectorType m_FirstLineElements, m_FirstColumnElements;
    RealVectorType m_FirstDiagonalElements, m_DiagonalElements, m_LastDiagonalElements;
    RealVectorType m_FirstLeftElements, m_FirstRightElements;
    double m_SecondLeftElement, m_SecondRightElement;

    RealVectorType m_FirstLineDerivativeElements, m_FirstColumnDerivativeElements;
    RealVectorType m_FirstDiagonalDerivativeElements, m_DiagonalDerivativeElements, m_LastDiagonalDerivativeElements;
    RealVectorType m_FirstLeftDerivativeElements, m_FirstRightDerivativeElements;
    double m_SecondLeftDerivativeElement, m_SecondRightDerivativeElement;

    // Internal work variables. Because of this, not thread safe !
    vnl_matrix <double> m_SimulatedT2Values;
    vnl_matrix <double> m_SimulatedDerivativeT2Values;
    RealVectorType m_OutputVector;
    RealVectorType m_OutputB1Derivative;
};
    
} // end namespace of anima
