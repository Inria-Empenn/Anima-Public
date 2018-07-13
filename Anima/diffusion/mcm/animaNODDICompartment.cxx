#include <animaNODDICompartment.h>
#include <animaVectorOperations.h>
#include <animaErrorFunctions.h>
#include <animaLevenbergTools.h>
#include <animaKummerFunctions.h>
#include <animaWatsonDistribution.h>
#include <animaMCMConstants.h>
#include <boost/math/special_functions/legendre.hpp>

namespace anima
{

void NODDICompartment::UpdateSignals(double bValue, const Vector3DType &gradient)
{
    if (std::abs(bValue - m_CurrentBValue) < 1.0e-6 && anima::ComputeNorm(gradient - m_CurrentGradient) < 1.0e-6 && !m_ModifiedParameters)
        return;
    
    // This implementation of NODDI signal expression follows three papers:
    // 1) Zhang et al., 2012, Neuroimage: paper that introduces NODDI, with
    // integral representation of the signal in terms of NODDI parameters
    // (Eqs 1-8).
    // 2) Zhang et al., 2011, Neuroimage: it shows how the integral
    // representation of the NODDI intra-axonal signal (Eq. A2 with L_\perp
    // = 0 and L_// = d) can be written in series expansion using the
    // spherical harmonic expansion of the Watson probability density (Eq.
    // A4).
    // 3) Jespersen et al., 2007, Neuroimage: it gives the analytic
    // expression of the C function involved in the SH series expansion of
    // NODDI intra-axonal signal (Eq. 14).
    
    this->UpdateKappaValues();
    
    double theta = this->GetOrientationTheta();
    double phi = this->GetOrientationPhi();
    double kappa = this->GetOrientationConcentration();
    double nuic = 1.0 - this->GetExtraAxonalFraction();
    double dpara = this->GetAxialDiffusivity();
    
    Vector3DType compartmentOrientation(0.0);
    anima::TransformSphericalToCartesianCoordinates(theta,phi,1.0,compartmentOrientation);
    double innerProd = anima::ComputeScalarProduct(gradient, compartmentOrientation);
    
    // Aic and its derivatives w.r.t. params
    m_IntraAxonalSignal = 0;
    m_IntraAngleDerivative = 0;
    m_IntraKappaDerivative = 0;
    m_IntraAxialDerivative = 0;
    double x = bValue * dpara;
    
    for (unsigned int i = 0;i < m_WatsonSHCoefficients.size();++i)
    {
        double coefVal = m_WatsonSHCoefficients[i];
        double sqrtVal = std::sqrt((4.0 * i + 1.0) / (4.0 * M_PI));
        double legendreVal = boost::math::legendre_p(2 * i, innerProd);
        double kummerVal = anima::KummerFunction(-x, i + 0.5, 2.0 * i + 1.5, false, true);
        double xPowVal = std::pow(-x, (double)i);
        double cVal = xPowVal * kummerVal;
        
        // Signal
        m_IntraAxonalSignal += coefVal * sqrtVal * legendreVal * cVal;
        
        // Derivatives
        double coefDerivVal = m_WatsonSHCoefficientDerivatives[i];
        double legendreDerivVal = boost::math::legendre_p_prime(2 * i, innerProd);
        
        double cDerivVal = 0.0;
        if (m_EstimateAxialDiffusivity)
        {
            cDerivVal = -xPowVal * anima::KummerFunction(-x, i + 1.5, 2.0 * i + 2.5, false, true);
            if (i > 0)
                cDerivVal += xPowVal * i * kummerVal / x;
        }
        
        m_IntraAngleDerivative += coefVal * sqrtVal * legendreDerivVal * cVal;
        m_IntraKappaDerivative += coefDerivVal * sqrtVal * legendreVal * cVal;
        m_IntraAxialDerivative += coefVal * sqrtVal * legendreVal * cDerivVal;
    }
    
    m_IntraAxonalSignal /= 2.0;
    m_IntraAngleDerivative /= 2.0;
    m_IntraKappaDerivative /= 2.0;
    m_IntraAxialDerivative /= 2.0;
    
    // Aec
    double appAxialDiff = dpara * (1.0 - nuic * (1.0 - m_Tau1));
    double appRadialDiff = dpara * (1.0 - nuic * (1.0 + m_Tau1) / 2.0);
    
    m_ExtraAxonalSignal = std::exp(-bValue * (appRadialDiff + (appAxialDiff - appRadialDiff) * innerProd * innerProd));
    
    m_CurrentBValue = bValue;
    m_CurrentGradient = gradient;
    m_ModifiedParameters = false;
}

double NODDICompartment::GetFourierTransformedDiffusionProfile(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient)
{
    double bValue = anima::GetBValueFromAcquisitionParameters(smallDelta, bigDelta, gradientStrength);
    this->UpdateSignals(bValue, gradient);
    double nuec = this->GetExtraAxonalFraction();
    double signal = (1.0 - nuec) * m_IntraAxonalSignal + nuec * m_ExtraAxonalSignal;
    return signal;
}
    
NODDICompartment::ListType &NODDICompartment::GetSignalAttenuationJacobian(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient)
{
    double bValue = anima::GetBValueFromAcquisitionParameters(smallDelta, bigDelta, gradientStrength);
    this->UpdateSignals(bValue, gradient);
    
    m_JacobianVector.resize(this->GetNumberOfParameters());
    
    double theta = this->GetOrientationTheta();
    double phi = this->GetOrientationPhi();
    double kappa = this->GetOrientationConcentration();
    double nuic = 1.0 - this->GetExtraAxonalFraction();
    double dpara = this->GetAxialDiffusivity();
    
    Vector3DType compartmentOrientation(0.0);
    anima::TransformSphericalToCartesianCoordinates(theta,phi,1.0,compartmentOrientation);
    double innerProd = anima::ComputeScalarProduct(gradient, compartmentOrientation);
    
    // Base extra angle derivative
    double extraAngleDerivative = -bValue * dpara * nuic * (3.0 * m_Tau1 - 1.0) * innerProd * m_ExtraAxonalSignal;
    
    //------------------------
    // Derivative w.r.t. theta
    //------------------------
    double thetaDeriv = 1.0;
    if (this->GetUseBoundedOptimization())
        thetaDeriv = levenberg::BoundedDerivativeAddOn(theta, this->GetBoundedSignVectorValue(0),
                                                       m_ZeroLowerBound, m_PolarAngleUpperBound);
    
    double innerProdThetaDeriv = (gradient[0] * std::cos(phi) + gradient[1] * std::sin(phi)) * std::cos(theta) - gradient[2] * std::sin(theta);
    
    // Derivative of Aic w.r.t. theta
    double intraThetaDerivative = m_IntraAngleDerivative * innerProdThetaDeriv;
    
    // Derivative of Aec w.r.t. theta
    double extraThetaDerivative = extraAngleDerivative * innerProdThetaDeriv;
    
    m_JacobianVector[0] = (nuic * intraThetaDerivative + (1.0 - nuic) * extraThetaDerivative) * thetaDeriv;
    
    //----------------------
    // Derivative w.r.t. phi
    //----------------------
    double phiDeriv = 1.0;
    if (this->GetUseBoundedOptimization())
        phiDeriv = levenberg::BoundedDerivativeAddOn(phi, this->GetBoundedSignVectorValue(1),
                                                     m_ZeroLowerBound, m_AzimuthAngleUpperBound);
    
    double innerProdPhiDeriv = std::sin(theta) * (gradient[1] * std::cos(phi) - gradient[0] * std::sin(phi));
    
    // Derivative of Aic w.r.t. phi
    double intraPhiDerivative = m_IntraAngleDerivative * innerProdPhiDeriv;
    
    // Derivative of Aec w.r.t. phi
    double extraPhiDerivative = extraAngleDerivative * innerProdPhiDeriv;
    
    m_JacobianVector[1] = (nuic * intraPhiDerivative + (1.0 - nuic) * extraPhiDerivative) * phiDeriv;
    
    //------------------------
    // Derivative w.r.t. kappa
    //------------------------
    unsigned int pos = 2;

    if (m_EstimateOrientationConcentration)
    {
        double kappaDeriv = 1.0;
        if (this->GetUseBoundedOptimization())
            kappaDeriv = levenberg::BoundedDerivativeAddOn(kappa, this->GetBoundedSignVectorValue(pos),
                                                           m_ZeroLowerBound, m_WatsonKappaUpperBound);

        // Derivative of Aec w.r.t. kappa
        double extraKappaDerivative = bValue * dpara * nuic * m_Tau1Deriv * (1.0 - 3.0 * innerProd * innerProd) * m_ExtraAxonalSignal / 2.0;

        m_JacobianVector[pos] = (nuic * m_IntraKappaDerivative + (1.0 - nuic) * extraKappaDerivative) * kappaDeriv;
        ++pos;
    }
    
    //---------------------
    // Derivative w.r.t. nu
    //---------------------
    double tmpVal = 1.0 + m_Tau1 - (3.0 * m_Tau1 - 1.0) * innerProd * innerProd;
    if (m_EstimateExtraAxonalFraction)
    {
        double nuDeriv = -1.0;
        if (this->GetUseBoundedOptimization())
            nuDeriv *= levenberg::BoundedDerivativeAddOn(1.0 - nuic, this->GetBoundedSignVectorValue(pos), m_ZeroLowerBound, m_FractionUpperBound);

        double extraNuicDeriv = bValue * dpara * tmpVal * m_ExtraAxonalSignal / 2.0;

        m_JacobianVector[pos] = (m_IntraAxonalSignal - m_ExtraAxonalSignal + (1.0 - nuic) * extraNuicDeriv) * nuDeriv;
        ++pos;
    }
    
    if (m_EstimateAxialDiffusivity)
    {
        //---------------------------
        // Derivative w.r.t. to dpara
        //---------------------------
        double dparaDeriv = 1.0;
        if (this->GetUseBoundedOptimization())
            dparaDeriv = levenberg::BoundedDerivativeAddOn(dpara,
                                                        this->GetBoundedSignVectorValue(pos),
                                                        m_ZeroLowerBound, m_DiffusivityUpperBound);
        
        m_JacobianVector[pos] = bValue * dparaDeriv * (nuic * m_IntraAxialDerivative - (1.0 - nuic) * m_ExtraAxonalSignal * (1.0 - nuic * tmpVal / 2.0));
    }
    
    return m_JacobianVector;
}

double NODDICompartment::GetLogDiffusionProfile(const Vector3DType &sample)
{
    // The PDF of the intra-axonal space is not well defined (degenerated from 3D to 1D).
    // So we use here the extra-axonal diffusion tensor to define a zero-mean 3D Gaussian
    // diffusion profile.
    
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

void NODDICompartment::SetOrientationTheta(double num)
{
    if (num != this->GetOrientationTheta())
    {
        m_ModifiedParameters = true;
        this->Superclass::SetOrientationTheta(num);
    }
}

void NODDICompartment::SetOrientationPhi(double num)
{
    if (num != this->GetOrientationPhi())
    {
        m_ModifiedParameters = true;
        this->Superclass::SetOrientationPhi(num);
    }
}

void NODDICompartment::SetOrientationConcentration(double num)
{
    if (num != this->GetOrientationConcentration())
    {
        m_ModifiedParameters = true;
        m_ModifiedConcentration = true;
        this->Superclass::SetOrientationConcentration(num);
    }
}

void NODDICompartment::SetExtraAxonalFraction(double num)
{
    if (num != this->GetExtraAxonalFraction())
    {
        m_ModifiedParameters = true;
        this->Superclass::SetExtraAxonalFraction(num);
    }
}

void NODDICompartment::SetAxialDiffusivity(double num)
{
    if (num != this->GetAxialDiffusivity())
    {
        m_ModifiedParameters = true;
        this->Superclass::SetAxialDiffusivity(num);
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

    unsigned int pos = 2;
    if (m_EstimateOrientationConcentration)
    {
        this->SetOrientationConcentration(m_BoundedVector[pos]);
        ++pos;
    }

    if (m_EstimateExtraAxonalFraction)
    {
        this->SetExtraAxonalFraction(m_BoundedVector[pos]);
        ++pos;
    }
    
    if (m_EstimateAxialDiffusivity)
        this->SetAxialDiffusivity(m_BoundedVector[pos]);
    else
    {
        // Constraint as in Zhang et al. 2012, Neuroimage.
        this->SetAxialDiffusivity(1.71e-3);
    }
}

NODDICompartment::ListType &NODDICompartment::GetParametersAsVector()
{
    m_ParametersVector.resize(this->GetNumberOfParameters());

    m_ParametersVector[0] = this->GetOrientationTheta();
    m_ParametersVector[1] = this->GetOrientationPhi();

    unsigned int pos = 2;
    if (m_EstimateOrientationConcentration)
    {
        m_ParametersVector[pos] = this->GetOrientationConcentration();
        ++pos;
    }

    if (m_EstimateExtraAxonalFraction)
    {
        m_ParametersVector[pos] = this->GetExtraAxonalFraction();
        ++pos;
    }
    
    if (m_EstimateAxialDiffusivity)
        m_ParametersVector[pos] = this->GetAxialDiffusivity();
    
    if (this->GetUseBoundedOptimization())
        this->UnboundParameters(m_ParametersVector);

    return m_ParametersVector;
}

NODDICompartment::ListType &NODDICompartment::GetParameterLowerBounds()
{
    m_ParametersLowerBoundsVector.resize(this->GetNumberOfParameters());
    
    std::fill(m_ParametersLowerBoundsVector.begin(),m_ParametersLowerBoundsVector.end(),m_ZeroLowerBound);
    
    if (m_EstimateAxialDiffusivity)
        m_ParametersLowerBoundsVector[this->GetNumberOfParameters() - 1] = m_DiffusivityLowerBound;
    
    return m_ParametersLowerBoundsVector;
}

NODDICompartment::ListType &NODDICompartment::GetParameterUpperBounds()
{
    m_ParametersUpperBoundsVector.resize(this->GetNumberOfParameters());

    m_ParametersUpperBoundsVector[0] = m_PolarAngleUpperBound;
    m_ParametersUpperBoundsVector[1] = m_AzimuthAngleUpperBound;

    unsigned int pos = 2;

    if (m_EstimateOrientationConcentration)
    {
        m_ParametersUpperBoundsVector[pos] = m_WatsonKappaUpperBound;
        ++pos;
    }

    if (m_EstimateExtraAxonalFraction)
    {
        m_ParametersUpperBoundsVector[pos] = m_FractionUpperBound;
        ++pos;
    }
    
    if (m_EstimateAxialDiffusivity)
        m_ParametersUpperBoundsVector[pos] = m_DiffusivityUpperBound;

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

    unsigned int pos = 2;
    if (m_EstimateOrientationConcentration)
    {
        m_BoundedVector[pos] = levenberg::ComputeBoundedValue(params[pos], inputSign, m_ZeroLowerBound, m_WatsonKappaUpperBound);
        this->SetBoundedSignVectorValue(pos,inputSign);
        ++pos;
    }

    if (m_EstimateExtraAxonalFraction)
    {
        m_BoundedVector[pos] = levenberg::ComputeBoundedValue(params[pos], inputSign, m_ZeroLowerBound, m_FractionUpperBound);
        this->SetBoundedSignVectorValue(pos,inputSign);
        ++pos;
    }
    
    if (m_EstimateAxialDiffusivity)
    {
        m_BoundedVector[pos] = levenberg::ComputeBoundedValue(params[pos], inputSign, m_DiffusivityLowerBound, m_DiffusivityUpperBound);
        this->SetBoundedSignVectorValue(pos,inputSign);
    }
}

void NODDICompartment::UnboundParameters(ListType &params)
{
    params[0] = levenberg::UnboundValue(params[0], m_ZeroLowerBound, m_PolarAngleUpperBound);
    params[1] = levenberg::UnboundValue(params[1], m_ZeroLowerBound, m_AzimuthAngleUpperBound);

    unsigned int pos = 2;
    if (m_EstimateOrientationConcentration)
    {
        params[pos] = levenberg::UnboundValue(params[pos], m_ZeroLowerBound, m_WatsonKappaUpperBound);
        ++pos;
    }

    if (m_EstimateExtraAxonalFraction)
    {
        params[pos] = levenberg::UnboundValue(params[pos], m_ZeroLowerBound, m_FractionUpperBound);
        ++pos;
    }
    
    if (m_EstimateAxialDiffusivity)
        params[pos] = levenberg::UnboundValue(params[pos], m_DiffusivityLowerBound, m_DiffusivityUpperBound);
}

void NODDICompartment::SetEstimateOrientationConcentration(bool arg)
{
    if (m_EstimateOrientationConcentration == arg)
        return;

    m_EstimateOrientationConcentration = arg;
    m_ChangedConstraints = true;
}

void NODDICompartment::SetEstimateAxialDiffusivity(bool arg)
{
    if (m_EstimateAxialDiffusivity == arg)
        return;

    m_EstimateAxialDiffusivity = arg;
    m_ChangedConstraints = true;
}

void NODDICompartment::SetEstimateExtraAxonalFraction(bool arg)
{
    if (m_EstimateExtraAxonalFraction == arg)
        return;

    m_EstimateExtraAxonalFraction = arg;
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
    
    // Display OD map (Zhang et al., 2012, Neuroimage, Eq. 9) to stay consistent
    // with the authors' original work but Would be better to show (3 tau1 - 1) / 2.
    double odi = compartmentVector[currentPos];
    this->SetOrientationConcentration(1.0 / std::tan(M_PI / 2.0 * odi));
    ++currentPos;
    
    // Display intra-axonal fraction map to stay consistent with NODDI
    // parametrization as in Zhang et al., 2012, Neuroimage.
    this->SetExtraAxonalFraction(1.0 - compartmentVector[currentPos]);
    ++currentPos;
    
    this->SetAxialDiffusivity(compartmentVector[currentPos]);
    
    m_ModifiedParameters = false;
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
    
    m_NumberOfParameters = 2;

    if (m_EstimateOrientationConcentration)
        ++m_NumberOfParameters;

    if (m_EstimateAxialDiffusivity)
        ++m_NumberOfParameters;

    if (m_EstimateExtraAxonalFraction)
        ++m_NumberOfParameters;
    
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
    
    // Display OD map (Zhang et al., 2012, Neuroimage, Eq. 9) to stay consistent
    // with the authors' original work but Would be better to show (3 tau1 - 1) / 2.
    m_CompartmentVector[currentPos] = 2.0 / M_PI * std::atan(1.0 / this->GetOrientationConcentration());
    ++currentPos;
    
    // Display intra-axonal fraction map to stay consistent with NODDI
    // parametrization as in Zhang et al., 2012, Neuroimage.
    m_CompartmentVector[currentPos] = 1.0 - this->GetExtraAxonalFraction();
    ++currentPos;
    
    m_CompartmentVector[currentPos] = this->GetAxialDiffusivity();
    
    return m_CompartmentVector;
}

const NODDICompartment::Matrix3DType &NODDICompartment::GetDiffusionTensor()
{
    // Get an apparent diffusion tensor based on the extra-axonal space.
    
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
    double dawsonValue = anima::EvaluateDawsonIntegral(std::sqrt(kappa), true);
    m_Tau1 = (1.0 / dawsonValue - 1.0) / (2.0 * kappa);
    m_Tau1Deriv = (1.0 - (1.0 - dawsonValue * (2.0 * kappa - 1.0)) / (2.0 * dawsonValue * dawsonValue)) / (2.0 * kappa * kappa);
    anima::GetStandardWatsonSHCoefficients(kappa,m_WatsonSHCoefficients,m_WatsonSHCoefficientDerivatives);
    
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
