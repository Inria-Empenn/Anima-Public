#include <animaNODDICompartment.h>
#include <animaVectorOperations.h>
#include <animaErrorFunctions.h>
#include <animaLevenbergTools.h>
#include <animaKummerFunctions.h>
#include <animaWatsonDistribution.h>
#include <boost/math/special_functions/legendre.hpp>

namespace anima
{
void NODDICompartment::GetIESignals(double bValue, const Vector3DType &gradient)
{
    double theta = this->GetOrientationTheta();
    double phi = this->GetOrientationPhi();
    double kappa = this->GetOrientationConcentration();
    double nuic = 1.0 - this->GetExtraAxonalFraction();
    double dpara = this->GetAxialDiffusivity();
    
    this->UpdateKappaValues();
    Vector3DType compartmentOrientation(0.0);
    anima::TransformSphericalToCartesianCoordinates(theta,phi,1.0,compartmentOrientation);
    double innerProd = anima::ComputeScalarProduct(gradient, compartmentOrientation);
    
    // Deal first with intra-axonal signal
    m_IntraAxonalSignal = 0.0;
    
    for (unsigned int i = 0;i < m_WatsonSHCoefficients.size();++i)
    {
        m_IntraAxonalSignal += m_WatsonSHCoefficients[i] * std::sqrt((4.0 * i + 1.0) / (4.0 * M_PI)) * boost::math::legendre_p(2 * i, innerProd) * std::pow(-bValue * dpara, (double)i) * std::tgamma(i + 0.5) / std::tgamma(2.0 * i + 1.5) * anima::KummerFunction(-bValue * dpara, i + 0.5, 2.0 * i + 1.5);
    }
    
    m_IntraAxonalSignal /= 2.0;
    
    // Now deal with extra-axonal signal
    double appAxialDiff = dpara * (1.0 - nuic * (1.0 - m_Tau1));
    double appRadialDiff = dpara * (1.0 - nuic * (1.0 + m_Tau1) / 2.0);
    
    m_ExtraAxonalSignal = std::exp(-bValue * (appRadialDiff + (appAxialDiff - appRadialDiff) * innerProd * innerProd));
}

double NODDICompartment::GetFourierTransformedDiffusionProfile(double bValue, const Vector3DType &gradient)
{
    this->GetIESignals(bValue, gradient);
    double nuic = 1.0 - this->GetExtraAxonalFraction();
    double signal = nuic * m_IntraAxonalSignal + (1.0 - nuic) * m_ExtraAxonalSignal;
    return signal;
}
    
NODDICompartment::ListType &NODDICompartment::GetSignalAttenuationJacobian(double bValue, const Vector3DType &gradient)
{
    this->GetIESignals(bValue, gradient);
    
    m_JacobianVector.resize(this->GetNumberOfParameters());
    
    double theta = this->GetOrientationTheta();
    double phi = this->GetOrientationPhi();
    double kappa = this->GetOrientationConcentration();
    double nuic = 1.0 - this->GetExtraAxonalFraction();
    double dpara = this->GetAxialDiffusivity();
    
    Vector3DType compartmentOrientation(0.0);
    anima::TransformSphericalToCartesianCoordinates(theta,phi,1.0,compartmentOrientation);
    double innerProd = anima::ComputeScalarProduct(gradient, compartmentOrientation);
    
    // Derivative w.r.t. theta
    double thetaDeriv = 1.0;
    if (this->GetUseBoundedOptimization())
        thetaDeriv = levenberg::BoundedDerivativeAddOn(theta, this->GetBoundedSignVectorValue(0),
                                                       m_ZeroLowerBound, m_PolarAngleUpperBound);
    
    double innerProdThetaDeriv = (gradient[0] * std::cos(phi) + gradient[1] * std::sin(phi)) * std::cos(theta) - gradient[2] * std::sin(theta);
    m_JacobianVector[0] = (kappa * nuic * m_IntegralForThetaDerivative - bValue * dpara * nuic * (1.0 - nuic) * (3.0 * m_Tau1 - 1.0) * innerProdThetaDeriv * innerProd * m_ExtraAxonalSignal) * thetaDeriv;
    
    // Derivative w.r.t. phi
    double phiDeriv = 1.0;
    if (this->GetUseBoundedOptimization())
        phiDeriv = levenberg::BoundedDerivativeAddOn(phi, this->GetBoundedSignVectorValue(1),
                                                     m_ZeroLowerBound, m_AzimuthAngleUpperBound);
    
    double innerProdPhiDeriv = std::sin(theta) * (gradient[1] * std::cos(phi) - gradient[0] * std::sin(phi));
    m_JacobianVector[1] = (kappa * nuic * m_IntegralForPhiDerivative - bValue * dpara * nuic * (1.0 - nuic) * (3.0 * m_Tau1 - 1.0) * innerProdPhiDeriv * innerProd * m_ExtraAxonalSignal) * phiDeriv;
    
    // Derivative w.r.t. kappa
    double kappaDeriv = 1.0;
    if (this->GetUseBoundedOptimization())
        kappaDeriv = levenberg::BoundedDerivativeAddOn(kappa, this->GetBoundedSignVectorValue(2),
                                                       m_ZeroLowerBound, m_WatsonKappaUpperBound);
    
    double intraKappaDeriv = m_IntegralForKappaDerivative - m_KummerRatio * m_IntraAxonalSignal;
    double extraKappaDeriv = bValue * dpara * nuic * m_Tau1Deriv * (1.0 - 3.0 * innerProd * innerProd) * m_ExtraAxonalSignal / 2.0;
    m_JacobianVector[2] = (nuic * intraKappaDeriv + (1.0 - nuic) * extraKappaDeriv) * kappaDeriv;
    
    // Derivative w.r.t. nu
    double nuDeriv = -1.0;
    if (this->GetUseBoundedOptimization())
        nuDeriv *= levenberg::BoundedDerivativeAddOn(kappa, this->GetBoundedSignVectorValue(3),
                                                       m_EAFLowerBound, m_EAFUpperBound);
    
    double tmpVal = 1.0 + m_Tau1 - (3.0 * m_Tau1 - 1.0) * innerProd * innerProd;
    double extraNuicDeriv = bValue * dpara * tmpVal / 2.0 * m_ExtraAxonalSignal;
    m_JacobianVector[3] = (m_IntraAxonalSignal - m_ExtraAxonalSignal + (1.0 - nuic) * extraNuicDeriv) * nuDeriv;
    
    if (m_EstimateAxialDiffusivity)
    {
        // Derivative w.r.t. to dpara
        double dparaDeriv = 1.0;
        if (this->GetUseBoundedOptimization())
            dparaDeriv = levenberg::BoundedDerivativeAddOn(dpara,
                                                        this->GetBoundedSignVectorValue(4),
                                                        m_ZeroLowerBound, m_DiffusivityUpperBound);
        
        m_JacobianVector[4] = -bValue * dparaDeriv * (nuic * m_IntegralForDparaDerivative + (1.0 - nuic) * m_ExtraAxonalSignal * (1.0 - nuic * tmpVal / 2.0));
    }
    
    return m_JacobianVector;
}

double NODDICompartment::GetLogDiffusionProfile(const Vector3DType &sample)
{
    this->UpdateKappaValues();
    
    double intraFraction = 1.0 - this->GetExtraAxonalFraction();
    double axialDiff = this->GetAxialDiffusivity();
    double appAxialDiff = axialDiff * (1.0 - intraFraction * (1.0 - m_Tau1));
    double appRadialDiff = axialDiff * (1.0 - intraFraction * (1.0 + m_Tau1) / 2.0);
    
    Vector3DType compartmentOrientation(0.0);
    anima::TransformSphericalToCartesianCoordinates(this->GetOrientationTheta(),this->GetOrientationPhi(),1.0,compartmentOrientation);
    
    double innerProd = anima::ComputeScalarProduct(compartmentOrientation,sample);
    
    double resVal = - 1.5 * std::log(2.0 * M_PI) - 0.5 * std::log(appAxialDiff) - std::log(appRadialDiff);
    resVal -= (sample.squared_magnitude() - (1.0 - appRadialDiff / appAxialDiff) * innerProd * innerProd) / (2.0 * appRadialDiff);
    
    return resVal;
}

void NODDICompartment::SetOrientationConcentration(double num)
{
    if (num != this->GetOrientationConcentration())
    {
        m_ModifiedConcentration = true;
        this->Superclass::SetOrientationConcentration(num);
    }
}

void NODDICompartment::SetParametersFromVector(const ListType &params)
{
    if (params.size() != this->GetNumberOfParameters())
        return;
    
    if (this->GetUseBoundedOptimization())
    {
        if (params.size() != this->GetBoundedSignVector().size())
            this->GetBoundedSignVector().resize(params.size());

        this->BoundParameters(params);
    }
    else
        m_BoundedVector = params;

    this->SetOrientationTheta(m_BoundedVector[0]);
    this->SetOrientationPhi(m_BoundedVector[1]);
    this->SetOrientationConcentration(m_BoundedVector[2]);
    this->SetExtraAxonalFraction(m_BoundedVector[3]);
    
    if (m_EstimateAxialDiffusivity)
        this->SetAxialDiffusivity(m_BoundedVector[4]);
    else
    {
        // Constraint as in Zhang et al. 2012, Neuroimage.
        this->SetAxialDiffusivity(1.7e-3);
    }
}

NODDICompartment::ListType &NODDICompartment::GetParametersAsVector()
{
    m_ParametersVector.resize(this->GetNumberOfParameters());

    m_ParametersVector[0] = this->GetOrientationTheta();
    m_ParametersVector[1] = this->GetOrientationPhi();
    m_ParametersVector[2] = this->GetOrientationConcentration();
    m_ParametersVector[3] = this->GetExtraAxonalFraction();
    
    if (m_EstimateAxialDiffusivity)
        m_ParametersVector[4] = this->GetAxialDiffusivity();
    
    if (this->GetUseBoundedOptimization())
        this->UnboundParameters(m_ParametersVector);

    return m_ParametersVector;
}

NODDICompartment::ListType &NODDICompartment::GetParameterLowerBounds()
{
    m_ParametersLowerBoundsVector.resize(this->GetNumberOfParameters());
    
    std::fill(m_ParametersLowerBoundsVector.begin(),m_ParametersLowerBoundsVector.end(),m_ZeroLowerBound);
    
    m_ParametersLowerBoundsVector[3] = m_EAFLowerBound;
    
    if (m_EstimateAxialDiffusivity)
        m_ParametersLowerBoundsVector[4] = m_DiffusivityLowerBound;
    
    return m_ParametersLowerBoundsVector;
}

NODDICompartment::ListType &NODDICompartment::GetParameterUpperBounds()
{
    m_ParametersUpperBoundsVector.resize(this->GetNumberOfParameters());

    m_ParametersUpperBoundsVector[0] = m_PolarAngleUpperBound;
    m_ParametersUpperBoundsVector[1] = m_AzimuthAngleUpperBound;
    m_ParametersUpperBoundsVector[2] = m_WatsonKappaUpperBound;
    m_ParametersUpperBoundsVector[3] = m_EAFUpperBound;
    
    if (m_EstimateAxialDiffusivity)
        m_ParametersUpperBoundsVector[4] = m_DiffusivityUpperBound;

    return m_ParametersUpperBoundsVector;
}
    
void NODDICompartment::BoundParameters(const ListType &params)
{
    m_BoundedVector.resize(params.size());
    
    double inputSign = 1;
    m_BoundedVector[0] = levenberg::ComputeBoundedValue(params[0], inputSign, m_ZeroLowerBound, m_PolarAngleUpperBound);
    this->SetBoundedSignVectorValue(0,inputSign);
    m_BoundedVector[1] = levenberg::ComputeBoundedValue(params[1], inputSign, m_ZeroLowerBound, m_AzimuthAngleUpperBound);
    this->SetBoundedSignVectorValue(1,inputSign);
    m_BoundedVector[2] = levenberg::ComputeBoundedValue(params[2], inputSign, m_ZeroLowerBound, m_WatsonKappaUpperBound);
    this->SetBoundedSignVectorValue(2,inputSign);
    m_BoundedVector[3] = levenberg::ComputeBoundedValue(params[3], inputSign, m_EAFLowerBound, m_EAFUpperBound);
    this->SetBoundedSignVectorValue(3,inputSign);
    
    if (m_EstimateAxialDiffusivity)
    {
        m_BoundedVector[4] = levenberg::ComputeBoundedValue(params[4], inputSign, m_DiffusivityLowerBound, m_DiffusivityUpperBound);
        this->SetBoundedSignVectorValue(4,inputSign);
    }
}

void NODDICompartment::UnboundParameters(ListType &params)
{
    params[0] = levenberg::UnboundValue(params[0], m_ZeroLowerBound, m_PolarAngleUpperBound);
    params[1] = levenberg::UnboundValue(params[1], m_ZeroLowerBound, m_AzimuthAngleUpperBound);
    params[2] = levenberg::UnboundValue(params[2], m_ZeroLowerBound, m_WatsonKappaUpperBound);
    params[3] = levenberg::UnboundValue(params[3], m_EAFLowerBound, m_EAFUpperBound);
    
    if (m_EstimateAxialDiffusivity)
        params[4] = levenberg::UnboundValue(params[4], m_DiffusivityLowerBound, m_DiffusivityUpperBound);
}

void NODDICompartment::SetEstimateAxialDiffusivity(bool arg)
{
    if (m_EstimateAxialDiffusivity == arg)
        return;
    
    m_EstimateAxialDiffusivity = arg;
    m_ChangedConstraints = true;
}

void NODDICompartment::SetCompartmentVector(ModelOutputVectorType &compartmentVector)
{
    if (compartmentVector.GetSize() != this->GetCompartmentSize())
        itkExceptionMacro("The input vector size does not match the size of the compartment");

    Vector3DType compartmentOrientation, sphDir;
    
    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        compartmentOrientation[i] = compartmentVector[i];
    
    anima::TransformCartesianToSphericalCoordinates(compartmentOrientation,sphDir);
    this->SetOrientationTheta(sphDir[0]);
    this->SetOrientationPhi(sphDir[1]);
    
    unsigned int currentPos = m_SpaceDimension;
    
    double odi = compartmentVector[currentPos];
    this->SetOrientationConcentration(1.0 / std::tan(M_PI / 2.0 * odi));
    ++currentPos;
    
    this->SetExtraAxonalFraction(1.0 - compartmentVector[currentPos]);
    ++currentPos;
    
    this->SetAxialDiffusivity(compartmentVector[currentPos]);
    
    m_ModifiedConcentration = false;
}

unsigned int NODDICompartment::GetCompartmentSize()
{
    return 6;
}

unsigned int NODDICompartment::GetNumberOfParameters()
{
    if (!m_ChangedConstraints)
        return m_NumberOfParameters;
    
    // The number of parameters before constraints is the compartment size - 2 because the orientation is stored in cartesian coordinates
    m_NumberOfParameters = this->GetCompartmentSize() - 1;
    
    if (!m_EstimateAxialDiffusivity)
        --m_NumberOfParameters;
    
    m_ChangedConstraints = false;
    return m_NumberOfParameters;
}

NODDICompartment::ModelOutputVectorType &NODDICompartment::GetCompartmentVector()
{
    if (m_CompartmentVector.GetSize() != this->GetCompartmentSize())
        m_CompartmentVector.SetSize(this->GetCompartmentSize());

    Vector3DType compartmentOrientation(0.0);
    anima::TransformSphericalToCartesianCoordinates(this->GetOrientationTheta(),this->GetOrientationPhi(),1.0,compartmentOrientation);
    
    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        m_CompartmentVector[i] = compartmentOrientation[i];
    
    unsigned int currentPos = m_SpaceDimension;
    
    m_CompartmentVector[currentPos] = 2.0 / M_PI * std::atan(1.0 / this->GetOrientationConcentration());
    ++currentPos;
    
    m_CompartmentVector[currentPos] = 1.0 - this->GetExtraAxonalFraction();
    ++currentPos;
    
    m_CompartmentVector[currentPos] = this->GetAxialDiffusivity();
    
    return m_CompartmentVector;
}

const NODDICompartment::Matrix3DType &NODDICompartment::GetDiffusionTensor()
{
    this->UpdateKappaValues();
    
    double intraFraction = 1.0 - this->GetExtraAxonalFraction();
    double axialDiff = this->GetAxialDiffusivity();
    double appAxialDiff = axialDiff * (1.0 - intraFraction * (1.0 - m_Tau1));
    double appRadialDiff = axialDiff * (1.0 - intraFraction * (1.0 + m_Tau1) / 2.0);
    
    m_DiffusionTensor.Fill(0.0);
    
    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        m_DiffusionTensor(i,i) = appRadialDiff;
    
    Vector3DType compartmentOrientation(0.0);
    anima::TransformSphericalToCartesianCoordinates(this->GetOrientationTheta(),this->GetOrientationPhi(),1.0,compartmentOrientation);
    
    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        for (unsigned int j = i;j < m_SpaceDimension;++j)
        {
            m_DiffusionTensor(i,j) += compartmentOrientation[i] * compartmentOrientation[j] * (appAxialDiff - appRadialDiff);
            if (i != j)
                m_DiffusionTensor(j,i) = m_DiffusionTensor(i,j);
        }
    
    return m_DiffusionTensor;
}

void NODDICompartment::UpdateKappaValues()
{
    if (!m_ModifiedConcentration)
        return;
    
    double kappa = this->GetOrientationConcentration();
    double kappaSqrt = std::sqrt(kappa);
    double kappaSq = kappa * kappa;
    
    double fVal = anima::EvaluateDawsonIntegral(kappaSqrt);
    m_Tau1 = -1.0 / (2.0 * kappa) + 1.0 / (2.0 * fVal * kappaSqrt);
    m_Tau1Deriv = 1.0 / (2.0 * kappaSq) - (1.0 - kappaSqrt * fVal * (2.0 - 1.0 / kappaSq)) / (4.0 * kappa * fVal * fVal);
    double kummerVal = anima::KummerFunction(kappa, 0.5, 1.5);
    m_KummerRatio = std::exp(kappa) * (1.0 - std::exp(-kappa) * kummerVal) / (2.0 * kappa * kummerVal);
    
    anima::GetStandardWatsonSHCoefficients(kappa,m_WatsonSHCoefficients);
    
    m_ModifiedConcentration = false;
}

double NODDICompartment::GetFractionalAnisotropy()
{
    double intraFraction = 1.0 - this->GetExtraAxonalFraction();
    double axialDiff = this->GetAxialDiffusivity();
    double l1 = axialDiff * (1.0 - intraFraction * (1.0 - m_Tau1));
    double l2 = axialDiff * (1.0 - intraFraction * (1.0 + m_Tau1) / 2.0);
    double l3 = l2;
    
    double numFA = std::sqrt ((l1 - l2) * (l1 - l2) + (l2 - l3) * (l2 - l3) + (l3 - l1) * (l3 - l1));
    double denomFA = std::sqrt (l1 * l1 + l2 * l2 + l3 * l3);
    
    double fa = 0;
    if (denomFA != 0)
        fa = std::sqrt(0.5) * (numFA / denomFA);
    
    return fa;
}

} //end namespace anima
