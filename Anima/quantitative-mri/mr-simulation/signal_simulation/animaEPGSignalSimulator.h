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

private:
    double m_EchoSpacing;
    double m_ExcitationFlipAngle;
    unsigned int m_NumberOfEchoes;

    double m_FirstEPGProduct, m_SecondEPGProduct, m_ThirdEPGProduct, m_FourthEPGProduct, m_FifthEPGProduct;
    double m_FirstDerivativeProduct, m_SecondDerivativeProduct, m_ThirdDerivativeProduct;

    // Internal work variables. Because of this, not thread safe !
    vnl_matrix <double> m_SimulatedT2Values;
    vnl_matrix <double> m_SimulatedDerivativeT2Values;
    RealVectorType m_OutputVector;
    RealVectorType m_OutputB1Derivative;
};
    
} // end namespace of anima
