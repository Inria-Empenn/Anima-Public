#include <animaStaniszCompartment.h>
#include <animaLevenbergTools.h>

#include <animaVectorOperations.h>
#include <itkSymmetricEigenAnalysis.h>
#include <animaMCMConstants.h>

namespace anima
{

const double StaniszCompartment::m_StaniszAxialDiffusivityLowerBound = 1.5e-4;

void StaniszCompartment::UpdateSignals(double smallDelta, double largeDelta, double gradientStrength, const Vector3DType &gradient)
{
    if (std::abs(smallDelta - m_CurrentSmallDelta) < 1.0e-6 && std::abs(largeDelta - m_CurrentLargeDelta) && std::abs(gradientStrength - m_CurrentGradientStrength) < 1.0e-6 && anima::ComputeNorm(gradient - m_CurrentGradient) < 1.0e-6 && !m_ModifiedParameters)
        return;

    double alpha = anima::DiffusionGyromagneticRatio * smallDelta * gradientStrength;
    double bValue = anima::GetBValueFromAcquisitionParameters(smallDelta, largeDelta, gradientStrength);
    double tissueRadius = this->GetTissueRadius();
    double axialDiff = this->GetAxialDiffusivity();
    double alphaRs = alpha * tissueRadius;
    double piSquare = M_PI * M_PI;
    double alphaRsSquare = alphaRs * alphaRs;
    double deltaDiff = largeDelta - smallDelta / 3.0;
    double epsilon = std::numeric_limits<double>::epsilon();

    m_FirstSummation = 0.0;
    m_SecondSummation = 0.0;
    m_ThirdSummation = 0.0;
    m_FourthSummation = 0.0;
    for (unsigned int i = 0;i < m_MaximumNumberOfSumElements;++i)
    {
        double n = i + 1.0;
        double npiSquare = n * n * piSquare;
        double denomValue = alphaRsSquare - npiSquare;
        
        if (denomValue < 1.0e-6)
            continue;
        
        double internalTermCos = std::exp(-npiSquare * deltaDiff * axialDiff / (tissueRadius * tissueRadius));
        double internalTermSin = internalTermCos;

        if ((i+1) % 2 == 0)
        {
            internalTermCos *= 1.0 - std::cos(alphaRs);
            internalTermSin *= std::sin(alphaRs);
        }
        else
        {
            internalTermCos *= 1.0 + std::cos(alphaRs);
            internalTermSin *= -std::sin(alphaRs);
        }

        internalTermCos /= denomValue * denomValue;
        internalTermSin /= denomValue * denomValue;

        m_FirstSummation += internalTermCos;
        m_SecondSummation += internalTermCos * n * n;
        m_ThirdSummation += internalTermSin;
        m_FourthSummation += internalTermCos / denomValue;

        bool stopSum1 = internalTermCos < epsilon * m_FirstSummation;
        bool stopSum2 = internalTermCos * n * n < epsilon * m_SecondSummation;
        bool stopSum3 = internalTermSin < epsilon * m_ThirdSummation;
        bool stopSum4 = internalTermCos < epsilon * m_FourthSummation * denomValue;

        if (stopSum1 && stopSum2 && stopSum3 && stopSum4)
            break;
    }

    m_CurrentSmallDelta = smallDelta;
    m_CurrentLargeDelta = largeDelta;
    m_CurrentGradientStrength = gradientStrength;
    m_CurrentGradient = gradient;
    m_ModifiedParameters = false;
}

double StaniszCompartment::GetFourierTransformedDiffusionProfile(double smallDelta, double largeDelta, double gradientStrength, const Vector3DType &gradient)
{
    this->UpdateSignals(smallDelta, largeDelta, gradientStrength, gradient);

    double alpha = anima::DiffusionGyromagneticRatio * smallDelta * gradientStrength;
    double alphaRs = alpha * this->GetTissueRadius();
    double alphaRsSquare = alphaRs * alphaRs;

    double signalValue = (alphaRs < 1.0e-6) ? 1.0 - alphaRsSquare / 12.0 : 2.0 * (1.0 - std::cos(alphaRs)) / (alphaRsSquare);
    signalValue += 4.0 * alphaRsSquare *  m_FirstSummation;

    return signalValue;
}

StaniszCompartment::ListType &StaniszCompartment::GetSignalAttenuationJacobian(double smallDelta, double largeDelta, double gradientStrength, const Vector3DType &gradient)
{
    this->UpdateSignals(smallDelta, largeDelta, gradientStrength, gradient);

    m_JacobianVector.resize(this->GetNumberOfParameters());

    double tissueRadius = this->GetTissueRadius();
    double axialDiff = this->GetAxialDiffusivity();
    double alpha = anima::DiffusionGyromagneticRatio * smallDelta * gradientStrength;
    double alphaRs = alpha * tissueRadius;
    double alphaRsSquare = alphaRs * alphaRs;
    double piSquare = M_PI * M_PI;
    double bValue = anima::GetBValueFromAcquisitionParameters(smallDelta, largeDelta, gradientStrength);

    // Derivative w.r.t. R_S
    m_JacobianVector[0] = (alphaRs < 1.0e-6) ? -alphaRsSquare / (6.0 * tissueRadius) : 2.0 * (alphaRs * std::sin(alphaRs) - 2.0 * (1.0 - std::cos(alphaRs))) / (alphaRsSquare * tissueRadius);
    m_JacobianVector[0] += 8.0 * alpha * alphaRs * m_FirstSummation;
    m_JacobianVector[0] += 8.0 * piSquare * bValue * axialDiff * m_SecondSummation / tissueRadius;
    m_JacobianVector[0] += 4.0 * alpha * alphaRsSquare * m_ThirdSummation;
    m_JacobianVector[0] -= 16.0 * alpha * alphaRs * alphaRsSquare * m_FourthSummation;

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

        m_JacobianVector[1] = - 4.0 * bValue * piSquare * m_SecondSummation * aDiffDeriv;
    }

    return m_JacobianVector;
}

double StaniszCompartment::GetLogDiffusionProfile(const Vector3DType &sample)
{
    throw itk::ExceptionObject(__FILE__, __LINE__, "Stanisz PDF not implemented yet", ITK_LOCATION);
}

void StaniszCompartment::SetTissueRadius(double num)
{
    if (num != this->GetTissueRadius())
    {
        m_ModifiedParameters = true;
        this->Superclass::SetTissueRadius(num);
    }
}

void StaniszCompartment::SetAxialDiffusivity(double num)
{
    if (num != this->GetAxialDiffusivity())
    {
        m_ModifiedParameters = true;
        this->Superclass::SetAxialDiffusivity(num);
    }
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
        this->SetAxialDiffusivity(m_BoundedVector[1]);
}

StaniszCompartment::ListType &StaniszCompartment::GetParametersAsVector()
{
    m_ParametersVector.resize(this->GetNumberOfParameters());

    m_ParametersVector[0] = this->GetTissueRadius();

    if (m_EstimateAxialDiffusivity)
        m_ParametersVector[1] = this->GetAxialDiffusivity();
    
    if (this->GetUseBoundedOptimization())
        this->UnboundParameters(m_ParametersVector);

    return m_ParametersVector;
}

StaniszCompartment::ListType &StaniszCompartment::GetParameterLowerBounds()
{
    m_ParametersLowerBoundsVector.resize(this->GetNumberOfParameters());

    m_ParametersLowerBoundsVector[0] = m_TissueRadiusLowerBound;

    if (m_EstimateAxialDiffusivity)
        m_ParametersLowerBoundsVector[1] = m_StaniszAxialDiffusivityLowerBound;

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
        m_BoundedVector[1] = levenberg::ComputeBoundedValue(params[1],inputSign, m_StaniszAxialDiffusivityLowerBound, m_DiffusivityUpperBound);
        this->SetBoundedSignVectorValue(1,inputSign);
    }
}

void StaniszCompartment::UnboundParameters(ListType &params)
{
    params[0] = levenberg::UnboundValue(params[0], m_TissueRadiusLowerBound, m_TissueRadiusUpperBound);
    
    if (m_EstimateAxialDiffusivity)
        params[1] = levenberg::UnboundValue(params[1], m_StaniszAxialDiffusivityLowerBound, m_DiffusivityUpperBound);
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

    m_ModifiedParameters = false;
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
