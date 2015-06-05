#pragma once

#include <complex>
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
    typedef std::vector < std::complex <double> > ComplexVectorType;
    typedef std::vector <ComplexVectorType> MatrixType;

    RealVectorType GetValue(double t1Value, double t2Value,
                            double b1Value, double m0Value);

    void SetEchoSpacing(double val) {m_EchoSpacing = val;}
    void SetExcitationFlipAngle(double val) {m_ExcitationFlipAngle = val;}
    void SetFlipAngle(double val) {m_FlipAngle = val;}

    void SetNumberOfEchoes(unsigned int val) {m_NumberOfEchoes = val;}

    void SetB1OnExcitationAngle(bool val) {m_B1OnExcitationAngle = val;}

protected:
    void ComputeT2SignalMatrices(double t1Value, double t2Value, double b1Value);

private:
    double m_EchoSpacing;
    double m_ExcitationFlipAngle;
    double m_FlipAngle;
    unsigned int m_NumberOfEchoes;

    bool m_B1OnExcitationAngle;

    std::vector <MatrixType> m_DiagonalT2SignalMatrices;
    MatrixType m_LowerDiagonalT2SignalMatrices;
    MatrixType m_UpperDiagonalT2SignalMatrices;
};
    
} // end namespace of anima
