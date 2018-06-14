#include <animaStaniszCompartment.h>
#include <animaLevenbergTools.h>

#include <animaVectorOperations.h>
#include <itkSymmetricEigenAnalysis.h>
#include <animaMCMConstants.h>

namespace anima
{

double StaniszCompartment::GetFourierTransformedDiffusionProfile(double smallDelta, double largeDelta, double gradientStrength, const Vector3DType &gradient)
{
    double alpha = anima::DiffusionGyromagneticRatio * smallDelta * gradientStrength;
    double bValue = anima::GetBValueFromAcquisitionParameters(smallDelta, largeDelta, gradientStrength);
    double tissueRadius = this->GetTissueRadius();
    double axialDiff = this->GetAxialDiffusivity();
    double alphaRs = alpha * tissueRadius;
    double piSquare = M_PI * M_PI;
    double alphaRsSquare = alphaRs * alphaRs;
    double deltaDiff = largeDelta - smallDelta / 3.0;

    double signalValue = (alphaRs < 1.0e-6) ? 1.0 - alphaRsSquare / 12.0 : 2.0 * (1.0 - std::cos(alphaRs)) / (alphaRsSquare);
    
    double sumValue = 0.0;
    for (unsigned int i = 0;i < m_MaximumNumberOfSumElements;++i)
    {
        double n = i + 1.0;
        double npiSquare = n * n * piSquare;
        double denomValue = alphaRsSquare - npiSquare;
        
        if (std::abs(denomValue) < 1.0e-6)
            continue;
        
        double internalTerm = std::exp(-npiSquare * deltaDiff * axialDiff / (tissueRadius * tissueRadius));

        if ((i+1) % 2 == 0)
            internalTerm *= 1.0 - std::cos(alphaRs);
        else
            internalTerm *= 1.0 + std::cos(alphaRs);

        internalTerm /= denomValue * denomValue;

        sumValue += internalTerm;

        if (internalTerm < std::numeric_limits<double>::epsilon() * sumValue)
            break;
    }

    signalValue += 4.0 * alphaRsSquare * sumValue;

    return signalValue;
}

StaniszCompartment::ListType &StaniszCompartment::GetSignalAttenuationJacobian(double smallDelta, double largeDelta, double gradientStrength, const Vector3DType &gradient)
{
    m_JacobianVector.resize(this->GetNumberOfParameters());

    double alpha = anima::DiffusionGyromagneticRatio * smallDelta * gradientStrength;
    double bValue = anima::GetBValueFromAcquisitionParameters(smallDelta, largeDelta, gradientStrength);
    double tissueRadius = this->GetTissueRadius();
    double axialDiff = this->GetAxialDiffusivity();
    double alphaRs = alpha * tissueRadius;
    double piSquare = M_PI * M_PI;
    double alphaRsSquare = alphaRs * alphaRs;
    double deltaDiff = largeDelta - smallDelta / 3.0;

    double derivativeSum2Value = 0.0;
    for (unsigned int i = 0;i < m_MaximumNumberOfSumElements;++i)
    {
        double n = i + 1.0;
        double npiSquare = n * n * piSquare;
        double denomValue = alphaRsSquare - npiSquare;
        
        if (std::abs(denomValue) < 1.0e-6)
            continue;
        
        double internalTerm = n * n * std::exp(-npiSquare * deltaDiff * axialDiff / (tissueRadius * tissueRadius));

        if ((i+1) % 2 == 0)
            internalTerm *= 1.0 - std::cos(alphaRs);
        else
            internalTerm *= 1.0 + std::cos(alphaRs);

        internalTerm /= denomValue * denomValue;

        derivativeSum2Value += internalTerm;

        if (internalTerm < std::numeric_limits<double>::epsilon() * derivativeSum2Value)
            break;
    }

    // Derivative w.r.t. R_S
    m_JacobianVector[0] = (alphaRs < 1.0e-6) ? -alphaRsSquare / (6.0 * tissueRadius) : 2.0 * (alphaRs * std::sin(alphaRs) - 2.0 * (1.0 - std::cos(alphaRs))) / (alphaRsSquare * tissueRadius);
    
    double derivativeSum1Value = 0.0;
    double derivativeSum3Value = 0.0;
    double derivativeSum4Value = 0.0;
    for (unsigned int i = 0;i < m_MaximumNumberOfSumElements;++i)
    {
        double n = i + 1.0;
        double npiSquare = n * n * piSquare;
        double alphaRsNPI = alphaRsSquare - npiSquare;
        
        if (std::abs(alphaRsNPI) < 1.0e-6)
            continue;
        
        double internalTermCos = std::exp(-npiSquare * deltaDiff * axialDiff / (tissueRadius * tissueRadius));
        double internalTermSine = internalTermCos;

        if ((i+1) % 2 == 0)
        {
            internalTermCos *= 1.0 - std::cos(alphaRs);
            internalTermSine *= std::sin(alphaRs);
        }
        else
        {
            internalTermCos *= 1.0 + std::cos(alphaRs);
            internalTermSine *= - std::sin(alphaRs);
        }

        internalTermCos /= alphaRsNPI * alphaRsNPI;
        internalTermSine /= alphaRsNPI * alphaRsNPI;

        derivativeSum1Value += internalTermCos;
        derivativeSum3Value += internalTermSine;
        derivativeSum4Value += internalTermCos / alphaRsNPI;

        if (internalTermCos < std::numeric_limits<double>::epsilon() * derivativeSum1Value)
            break;
        
        if (internalTermSine < std::numeric_limits<double>::epsilon() * derivativeSum3Value)
            break;
        
        if (internalTermCos < std::numeric_limits<double>::epsilon() * derivativeSum4Value * alphaRsNPI)
            break;
    }

    m_JacobianVector[0] += 8.0 * alpha * alphaRs * derivativeSum1Value;
    m_JacobianVector[0] += 8.0 * piSquare * bValue * axialDiff * derivativeSum2Value / tissueRadius;
    m_JacobianVector[0] += 4.0 * alpha * alphaRsSquare * derivativeSum3Value;
    m_JacobianVector[0] -= 16.0 * alpha * alphaRs * alphaRsSquare * derivativeSum4Value;

    double rsDeriv = 1.0;
    if (this->GetUseBoundedOptimization())
        rsDeriv = levenberg::BoundedDerivativeAddOn(tissueRadius, this->GetBoundedSignVectorValue(0),
                                                    m_TissueRadiusLowerBound, m_TissueRadiusUpperBound);

    m_JacobianVector[0] *= rsDeriv;

    if (m_EstimateAxialDiffusivity)
    {
    	// Derivative w.r.t. D_I
	    double aDiffDeriv = 1.0;
	    if (this->GetUseBoundedOptimization())
	        aDiffDeriv = levenberg::BoundedDerivativeAddOn(axialDiff, this->GetBoundedSignVectorValue(1),
	                                                       m_DiffusivityLowerBound, m_DiffusivityUpperBound);

	    m_JacobianVector[1] = - 4.0 * bValue * piSquare * derivativeSum2Value * aDiffDeriv;
    }

    return m_JacobianVector;
}

double StaniszCompartment::GetLogDiffusionProfile(const Vector3DType &sample)
{
    throw itk::ExceptionObject(__FILE__, __LINE__, "Stanisz PDF not implemented yet", ITK_LOCATION);
}

void StaniszCompartment::SetParametersFromVector(const ListType &params)
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

    this->SetTissueRadius(m_BoundedVector[0]);

    if (m_EstimateAxialDiffusivity)
        this->SetAxialDiffusivity(m_BoundedVector[1] + this->GetRadialDiffusivity1());
}

StaniszCompartment::ListType &StaniszCompartment::GetParametersAsVector()
{
    m_ParametersVector.resize(this->GetNumberOfParameters());

    m_ParametersVector[0] = this->GetTissueRadius();

    if (m_EstimateAxialDiffusivity)
    	m_ParametersVector[1] = this->GetAxialDiffusivity() - this->GetRadialDiffusivity1();
    
    if (this->GetUseBoundedOptimization())
        this->UnboundParameters(m_ParametersVector);

    return m_ParametersVector;
}

StaniszCompartment::ListType &StaniszCompartment::GetParameterLowerBounds()
{
    m_ParametersLowerBoundsVector.resize(this->GetNumberOfParameters());

    m_ParametersLowerBoundsVector[0] = m_TissueRadiusLowerBound;

    if (m_EstimateAxialDiffusivity)
    	m_ParametersLowerBoundsVector[1] = m_ZeroLowerBound;

    return m_ParametersLowerBoundsVector;
}

StaniszCompartment::ListType &StaniszCompartment::GetParameterUpperBounds()
{
    m_ParametersUpperBoundsVector.resize(this->GetNumberOfParameters());

    m_ParametersUpperBoundsVector[0] = m_TissueRadiusUpperBound;

    if (m_EstimateAxialDiffusivity)
    	m_ParametersUpperBoundsVector[1] = m_DiffusivityUpperBound;

    return m_ParametersUpperBoundsVector;
}

void StaniszCompartment::BoundParameters(const ListType &params)
{
    m_BoundedVector.resize(params.size());
    
    double inputSign = 1;
    m_BoundedVector[0] = levenberg::ComputeBoundedValue(params[0], inputSign, m_TissueRadiusLowerBound, m_TissueRadiusUpperBound);
    this->SetBoundedSignVectorValue(0,inputSign);
    
    if (m_EstimateAxialDiffusivity)
    {
        m_BoundedVector[1] = levenberg::ComputeBoundedValue(params[1],inputSign, m_ZeroLowerBound, m_DiffusivityUpperBound);
        this->SetBoundedSignVectorValue(1,inputSign);
    }
}

void StaniszCompartment::UnboundParameters(ListType &params)
{
    params[0] = levenberg::UnboundValue(params[0], m_TissueRadiusLowerBound, m_TissueRadiusUpperBound);
    
    if (m_EstimateAxialDiffusivity)
        params[1] = levenberg::UnboundValue(params[1], m_ZeroLowerBound, m_DiffusivityUpperBound);
}

void StaniszCompartment::SetEstimateAxialDiffusivity(bool arg)
{
    if (m_EstimateAxialDiffusivity == arg)
        return;

    m_EstimateAxialDiffusivity = arg;
    m_ChangedConstraints = true;
}

void StaniszCompartment::SetCompartmentVector(ModelOutputVectorType &compartmentVector)
{
    if (compartmentVector.GetSize() != this->GetCompartmentSize())
        itkExceptionMacro("The input vector size does not match the size of the compartment");

    this->SetTissueRadius(compartmentVector[0]);
    this->SetAxialDiffusivity(compartmentVector[1]);
}

unsigned int StaniszCompartment::GetCompartmentSize()
{
    return 2;
}

unsigned int StaniszCompartment::GetNumberOfParameters()
{
    if (!m_ChangedConstraints)
        return m_NumberOfParameters;

    // The number of parameters before constraints is 1 because we assume for now the radius is always estimated
    m_NumberOfParameters = 2;

    if (!m_EstimateAxialDiffusivity)
        --m_NumberOfParameters;

    m_ChangedConstraints = false;
    return m_NumberOfParameters;
}

StaniszCompartment::ModelOutputVectorType &StaniszCompartment::GetCompartmentVector()
{
    if (m_CompartmentVector.GetSize() != this->GetCompartmentSize())
        m_CompartmentVector.SetSize(this->GetCompartmentSize());

    m_CompartmentVector[0] = this->GetTissueRadius();
    m_CompartmentVector[1] = this->GetAxialDiffusivity();

    return m_CompartmentVector;
}

double StaniszCompartment::GetFractionalAnisotropy()
{
    return 0.0;
}

} //end namespace anima
