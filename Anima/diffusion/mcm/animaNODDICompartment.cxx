#include <animaNODDICompartment.h>
#include <animaVectorOperations.h>
#include <animaErrorFunctions.h>
#include <animaDistributionSampling.h>
#include <animaLevenbergTools.h>
#include <animaKummerFunctions.h>

namespace anima
{
double NODDICompartment::GetFourierTransformedDiffusionProfile(double bValue, const Vector3DType &gradient)
{
    this->UpdateIESignals(bValue, gradient);
    double nuic = 1.0 - this->GetExtraAxonalFraction();
    return nuic * m_IntraAxonalSignal + (1.0 - nuic) * m_ExtraAxonalSignal;
}
    
NODDICompartment::ListType &NODDICompartment::GetSignalAttenuationJacobian(double bValue, const Vector3DType &gradient)
{
    this->UpdateIESignals(bValue, gradient);
    
    m_JacobianVector.resize(this->GetNumberOfParameters());
    
    double kappa = this->GetOrientationConcentration();
    double nuic = 1.0 - this->GetExtraAxonalFraction();
    double dpara = this->GetAxialDiffusivity();
    double phi = this->GetOrientationPhi();
    double theta = this->GetOrientationTheta();
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
    
    double kummerNumerator = anima::KummerFunction(kappa, 1.5, 2.5);
    double kummerDenominator = anima::KummerFunction(kappa, 0.5, 1.5);
    double intraKappaDeriv = m_IntegralForKappaDerivative - kummerNumerator * m_IntraAxonalSignal / (3.0 * kummerDenominator);
    double extraKappaDeriv = bValue * dpara * nuic * m_Tau1Deriv * (1.0 - 3.0 * innerProd * innerProd) * m_ExtraAxonalSignal / 2.0;
    m_JacobianVector[2] = (nuic * intraKappaDeriv + (1.0 - nuic) * extraKappaDeriv) * kappaDeriv;
    
    // Derivative w.r.t. nu
    double nuDeriv = -1.0;
    if (this->GetUseBoundedOptimization())
        nuDeriv *= levenberg::BoundedDerivativeAddOn(kappa, this->GetBoundedSignVectorValue(3),
                                                       m_EAFLowerBound, m_EAFUpperBound);
    
    double tmpVal = 1.0 + m_Tau1 - (3.0 * m_Tau1 - 1.0) * innerProd * innerProd;
    double extraNuicDeriv = bValue * dpara * tmpVal / 2.0;
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
    return 1;
}

void NODDICompartment::SetOrientationTheta(double num)
{
    if (num != this->GetOrientationTheta())
    {
        m_ModifiedTheta = true;
        this->Superclass::SetOrientationTheta(num);
    }
}

void NODDICompartment::SetOrientationPhi(double num)
{
    if (num != this->GetOrientationPhi())
    {
        m_ModifiedPhi = true;
        this->Superclass::SetOrientationPhi(num);
    }
}

void NODDICompartment::SetOrientationConcentration(double num)
{
    if (num != this->GetOrientationConcentration())
    {
        m_ModifiedConcentration = true;
        this->Superclass::SetOrientationConcentration(num);
    }
}

void NODDICompartment::SetExtraAxonalFraction(double num)
{
    if (num != this->GetExtraAxonalFraction())
    {
        m_ModifiedEAF = true;
        this->Superclass::SetExtraAxonalFraction(num);
    }
}

void NODDICompartment::SetAxialDiffusivity(double num)
{
    if (num != this->GetAxialDiffusivity())
    {
        m_ModifiedAxialDiffusivity = true;
        this->Superclass::SetAxialDiffusivity(num);
    }
}

void NODDICompartment::SetRadialDiffusivity1(double num)
{
    this->Superclass::SetRadialDiffusivity1(num);
    this->SetRadialDiffusivity2(num);
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
    
    // Tortuosity model as in Zhang et al. 2012, Neuroimage
    this->SetRadialDiffusivity1(this->GetExtraAxonalFraction() * this->GetAxialDiffusivity());
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
    
    m_ModifiedTheta = false;
    m_ModifiedPhi = false;
    m_ModifiedConcentration = false;
    m_ModifiedEAF = false;
    m_ModifiedAxialDiffusivity = false;
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
    this->UpdateWatsonSamples();
    
    double axialDiff = this->GetAxialDiffusivity() * (1.0 - (1.0 - this->GetExtraAxonalFraction()) * (1.0 - m_Tau1));
    double radialDiff = this->GetAxialDiffusivity() * (1.0 - (1.0 - this->GetExtraAxonalFraction()) * (1.0 + m_Tau1) / 2.0);
    
    m_DiffusionTensor.Fill(0.0);
    
    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        m_DiffusionTensor(i,i) = radialDiff;
    
    Vector3DType compartmentOrientation(0.0);
    anima::TransformSphericalToCartesianCoordinates(this->GetOrientationTheta(),this->GetOrientationPhi(),1.0,compartmentOrientation);
    
    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        for (unsigned int j = i;j < m_SpaceDimension;++j)
        {
            m_DiffusionTensor(i,j) += compartmentOrientation[i] * compartmentOrientation[j] * (axialDiff - radialDiff);
            if (i != j)
                m_DiffusionTensor(j,i) = m_DiffusionTensor(i,j);
        }
    
    return m_DiffusionTensor;
}
    
void NODDICompartment::UpdateWatsonSamples()
{
    if (!m_ModifiedConcentration)
        return;
    
    double kappa = this->GetOrientationConcentration();
    double kappaSqrt = std::sqrt(kappa);
    double kappaSq = kappa * kappa;
    
    double fVal = anima::EvaluateDawsonFunction(kappaSqrt);
    m_Tau1 = -1.0 / (2.0 * kappa) + 1.0 / (2.0 * fVal * kappaSqrt);
    m_Tau1Deriv = 1.0 / (2.0 * kappaSq) - (1.0 - kappaSqrt * fVal * (2.0 - 1.0 / kappaSq)) / (4.0 * kappa * fVal * fVal);
    
    m_WatsonSamples.resize(m_NumberOfSamples);
    std::mt19937 generator(time(0));
    
    for (unsigned int i = 0;i < m_NumberOfSamples;++i)
        anima::SampleFromWatsonDistribution(kappa, m_NorthPole, m_WatsonSamples[i], 3, generator);
    
    m_ModifiedConcentration = false;
}

void NODDICompartment::UpdateIESignals(double bValue, const Vector3DType &gradient)
{
    if (!m_ModifiedTheta && !m_ModifiedPhi && !m_ModifiedConcentration && !m_ModifiedEAF && !m_ModifiedAxialDiffusivity)
        return;
    
    this->UpdateWatsonSamples();
    
    Vector3DType compartmentOrientation(0.0), tmpVec = m_NorthPole;
    
    anima::TransformSphericalToCartesianCoordinates(this->GetOrientationTheta(),this->GetOrientationPhi(),1.0,compartmentOrientation);
    
    Vector3DType rotationNormal;
    anima::ComputeCrossProduct(compartmentOrientation, m_NorthPole, rotationNormal);
    anima::Normalize(rotationNormal, rotationNormal);
    anima::RotateAroundAxis(gradient, this->GetOrientationTheta(), rotationNormal, tmpVec);
    
    double axialDiff = this->GetAxialDiffusivity();
    double radialDiff = this->GetRadialDiffusivity1();
    double appAxialDiff = axialDiff - (axialDiff - radialDiff) * (1.0 - m_Tau1);
    double appRadialDiff = axialDiff - (axialDiff - radialDiff) * (1.0 + m_Tau1) / 2.0;
    
    m_IntraAxonalSignal = 0;
    m_IntegralForDparaDerivative = 0;
    m_IntegralForKappaDerivative = 0;
    m_IntegralForThetaDerivative = 0;
    m_IntegralForPhiDerivative = 0;
    
    for (unsigned int i = 0;i < m_NumberOfSamples;++i)
    {
        double innerProd = anima::ComputeScalarProduct(m_WatsonSamples[i], tmpVec);
        double expVal = std::exp(-bValue * axialDiff * innerProd * innerProd);
        m_IntraAxonalSignal += expVal;
        m_IntegralForDparaDerivative += expVal * innerProd * innerProd;
        m_IntegralForKappaDerivative += expVal * m_WatsonSamples[i][2] * m_WatsonSamples[i][2];
        // Apply inverse rotation to Watson sample
        anima::RotateAroundAxis(m_WatsonSamples[i], -this->GetOrientationTheta(), rotationNormal, tmpVec);
        double angleVal = 2.0 * m_WatsonSamples[i][2];
        double thetaVal = angleVal * ((tmpVec[0] * std::cos(this->GetOrientationPhi()) + tmpVec[1] * std::sin(this->GetOrientationPhi())) * std::cos(this->GetOrientationTheta()) - tmpVec[2] * std::sin(this->GetOrientationTheta()));
        double phiVal = angleVal * std::sin(this->GetOrientationTheta()) * (tmpVec[1] * std::cos(this->GetOrientationPhi()) - tmpVec[0] * std::sin(this->GetOrientationPhi()));
        m_IntegralForThetaDerivative += expVal * thetaVal;
        m_IntegralForPhiDerivative += expVal * phiVal;
    }
    
    m_IntraAxonalSignal /= m_NumberOfSamples;
    m_IntegralForDparaDerivative /= m_NumberOfSamples;
    m_IntegralForKappaDerivative /= m_NumberOfSamples;
    m_IntegralForThetaDerivative /= m_NumberOfSamples;
    m_IntegralForPhiDerivative /= m_NumberOfSamples;
    
    double scalProd = anima::ComputeScalarProduct(gradient, compartmentOrientation);
    m_ExtraAxonalSignal = std::exp(-bValue * (appRadialDiff + (appAxialDiff - appRadialDiff) * scalProd * scalProd));
    
    m_ModifiedTheta = false;
    m_ModifiedPhi = false;
    m_ModifiedConcentration = false;
    m_ModifiedEAF = false;
    m_ModifiedAxialDiffusivity = false;
}

} //end namespace anima
