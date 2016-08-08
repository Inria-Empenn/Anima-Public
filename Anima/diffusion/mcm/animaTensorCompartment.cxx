#include <animaTensorCompartment.h>

#include <itkSymmetricEigenAnalysis.h>
#include <animaBaseTensorTools.h>

namespace anima
{

double TensorCompartment::GetFourierTransformedDiffusionProfile(double bValue, const Vector3DType &gradient)
{
    this->UpdateDiffusionTensor();

    double quadForm = 0;

    for (unsigned int i = 0;i < m_SpaceDimension;++i)
    {
        quadForm += m_DiffusionTensor(i,i) * gradient[i] * gradient[i];
        for (unsigned int j = i+1;j < m_SpaceDimension;++j)
            quadForm += 2 * m_DiffusionTensor(i,j) * gradient[i] * gradient[j];
    }

    return std::exp(- bValue * quadForm);
}
    
TensorCompartment::ListType TensorCompartment::GetSignalAttenuationJacobian(double bValue, const Vector3DType &gradient)
{
    this->UpdateDiffusionTensor();
    
    ListType jacobian(this->GetNumberOfParameters(),0.0);
    
    double signalAttenuation = this->GetFourierTransformedDiffusionProfile(bValue, gradient);
    double innerProd1 = anima::ComputeScalarProduct(gradient, m_EigenVector1);
    double innerProd2 = anima::ComputeScalarProduct(gradient, m_EigenVector2);
    
    double DgTe1DTheta = m_CosTheta * (gradient[0] * m_CosPhi + gradient[1] * m_SinPhi) - gradient[2] * m_SinTheta;
    double DgTe1DPhi = m_SinTheta * (gradient[1] * m_CosPhi - gradient[0] * m_SinPhi);
    
    double DgTe2DTheta = m_SinAlpha * innerProd1;
    double DgTe2DPhi = gradient[0] * (m_CosTheta * m_SinPhi * m_SinAlpha - m_CosPhi * m_CosAlpha) - gradient[1] * (m_SinPhi * m_CosAlpha + m_CosTheta * m_CosPhi * m_SinAlpha);
    double DgTe2DAlpha = gradient[0] * (m_SinPhi * m_SinAlpha - m_CosTheta * m_CosPhi * m_CosAlpha) - gradient[1] * (m_CosPhi * m_SinAlpha + m_CosTheta * m_SinPhi * m_CosAlpha) + gradient[2] * m_SinTheta * m_CosAlpha;
    
    // Derivative w.r.t. theta
    double thetaDeriv = 1.0;
    if (this->GetUseBoundedOptimization())
        thetaDeriv = mcm_utilities::BoundedDerivativeAddOn(this->GetOrientationTheta(), this->GetBoundedSignVectorValue(0),
                                                           m_ZeroLowerBound, m_PolarAngleUpperBound);

    jacobian[0] = -2.0 * bValue * ((this->GetAxialDiffusivity() - this->GetRadialDiffusivity2()) * innerProd1 * DgTe1DTheta + (this->GetRadialDiffusivity1() - this->GetRadialDiffusivity2()) * innerProd2 * DgTe2DTheta) * signalAttenuation * thetaDeriv;
    
    // Derivative w.r.t. phi
    double phiDeriv = 1.0;
    if (this->GetUseBoundedOptimization())
        phiDeriv = mcm_utilities::BoundedDerivativeAddOn(this->GetOrientationPhi(), this->GetBoundedSignVectorValue(1),
                                                         m_ZeroLowerBound, m_AzimuthAngleUpperBound);

    jacobian[1] = -2.0 * bValue * ((this->GetAxialDiffusivity() - this->GetRadialDiffusivity2()) * innerProd1 * DgTe1DPhi + (this->GetRadialDiffusivity1() - this->GetRadialDiffusivity2()) * innerProd2 * DgTe2DPhi) * signalAttenuation * phiDeriv;
    
    // Derivative w.r.t. alpha
    double alphaDeriv = 1.0;
    if (this->GetUseBoundedOptimization())
        alphaDeriv = mcm_utilities::BoundedDerivativeAddOn(this->GetPerpendicularAngle(), this->GetBoundedSignVectorValue(2),
                                                           m_ZeroLowerBound, m_AzimuthAngleUpperBound);

    jacobian[2] = -2.0 * bValue * (this->GetRadialDiffusivity1() - this->GetRadialDiffusivity2()) * innerProd2 * DgTe2DAlpha * signalAttenuation * alphaDeriv;

    if (m_EstimateDiffusivities)
    {
        // Derivative w.r.t. to d1
        double d1Deriv = 1.0;
        if (this->GetUseBoundedOptimization())
            d1Deriv = mcm_utilities::BoundedDerivativeAddOn(this->GetAxialDiffusivity() - this->GetRadialDiffusivity1(),
                                                            this->GetBoundedSignVectorValue(3),
                                                            m_ZeroLowerBound, m_DiffusivityUpperBound);
        
        jacobian[3] = -bValue * innerProd1 * innerProd1 * signalAttenuation * d1Deriv;
        
        // Derivative w.r.t. to d2
        double d2Deriv = 1.0;
        if (this->GetUseBoundedOptimization())
            d2Deriv = mcm_utilities::BoundedDerivativeAddOn(this->GetRadialDiffusivity1() - this->GetRadialDiffusivity2(),
                                                            this->GetBoundedSignVectorValue(4),
                                                            m_ZeroLowerBound, m_DiffusivityUpperBound);

        jacobian[4] = -bValue * (innerProd1 * innerProd1 + innerProd2 * innerProd2) * signalAttenuation * d2Deriv;
        
        // Derivative w.r.t. to d3
        double d3Deriv = 1.0;
        if (this->GetUseBoundedOptimization())
            d3Deriv = mcm_utilities::BoundedDerivativeAddOn(this->GetRadialDiffusivity2(),
                                                            this->GetBoundedSignVectorValue(5),
                                                            m_DiffusivityLowerBound, m_RadialDiffusivityUpperBound);
        
        jacobian[5] = -bValue * signalAttenuation * d3Deriv;
    }
    
    return jacobian;
}

double TensorCompartment::GetLogDiffusionProfile(const Vector3DType &sample)
{
    this->UpdateDiffusionTensor();

    double resVal = - 1.5 * std::log(2.0 * M_PI) - 0.5 * std::log(m_TensorDeterminant);

    double quadForm = 0;

    for (unsigned int i = 0;i < m_SpaceDimension;++i)
    {
        quadForm += m_InverseDiffusionTensor(i,i) * sample[i] * sample[i];
        for (unsigned int j = i+1;j < m_SpaceDimension;++j)
            quadForm += 2 * m_InverseDiffusionTensor(i,j) * sample[i] * sample[j];
    }

    resVal -= quadForm / 2.0;

    return resVal;
}

void TensorCompartment::SetOrientationTheta(double num)
{
    if (num != this->GetOrientationTheta())
    {
        m_ModifiedTensor = true;
        m_ModifiedAngles = true;
        this->Superclass::SetOrientationTheta(num);
    }
}

double TensorCompartment::GetOrientationTheta()
{
    this->UpdateParametersFromTensor();
    return Superclass::GetOrientationTheta();
}

void TensorCompartment::SetOrientationPhi(double num)
{
    if (num != this->GetOrientationPhi())
    {
        m_ModifiedTensor = true;
        m_ModifiedAngles = true;
        this->Superclass::SetOrientationPhi(num);
    }
}

double TensorCompartment::GetOrientationPhi()
{
    this->UpdateParametersFromTensor();
    return Superclass::GetOrientationPhi();
}

void TensorCompartment::SetPerpendicularAngle(double num)
{
    if (num != this->GetPerpendicularAngle())
    {
        m_ModifiedTensor = true;
        m_ModifiedAngles = true;
        this->Superclass::SetPerpendicularAngle(num);
    }
}

double TensorCompartment::GetPerpendicularAngle()
{
    this->UpdateParametersFromTensor();
    return Superclass::GetPerpendicularAngle();
}

void TensorCompartment::SetAxialDiffusivity(double num)
{
    if (num != this->GetAxialDiffusivity())
    {
        m_ModifiedTensor = true;
        this->Superclass::SetAxialDiffusivity(num);
    }
}

double TensorCompartment::GetAxialDiffusivity()
{
    this->UpdateParametersFromTensor();
    return Superclass::GetAxialDiffusivity();
}

void TensorCompartment::SetRadialDiffusivity1(double num)
{
    if (num != this->GetRadialDiffusivity1())
    {
        m_ModifiedTensor = true;
        this->Superclass::SetRadialDiffusivity1(num);
    }
}

double TensorCompartment::GetRadialDiffusivity1()
{
    this->UpdateParametersFromTensor();
    return Superclass::GetRadialDiffusivity1();
}

void TensorCompartment::SetRadialDiffusivity2(double num)
{
    if (num != this->GetRadialDiffusivity2())
    {
        m_ModifiedTensor = true;
        this->Superclass::SetRadialDiffusivity2(num);
    }
}

double TensorCompartment::GetRadialDiffusivity2()
{
    this->UpdateParametersFromTensor();
    return Superclass::GetRadialDiffusivity2();
}

void TensorCompartment::SetParametersFromVector(const ListType &params)
{
    if (params.size() != this->GetNumberOfParameters())
        return;
    
    ListType boundedParams;
    if (this->GetUseBoundedOptimization())
    {
        if (params.size() != this->GetBoundedSignVector().size())
            this->GetBoundedSignVector().resize(params.size());

        boundedParams = this->BoundParameters(params);
    }
    else
        boundedParams = params;

    this->SetOrientationTheta(boundedParams[0]);
    this->SetOrientationPhi(boundedParams[1]);
    this->SetPerpendicularAngle(boundedParams[2]);

    if (m_EstimateDiffusivities)
    {
        this->SetAxialDiffusivity(boundedParams[3] + boundedParams[4] + boundedParams[5]);
        this->SetRadialDiffusivity1(boundedParams[4] + boundedParams[5]);
        this->SetRadialDiffusivity2(boundedParams[5]);
    }
}

TensorCompartment::ListType TensorCompartment::GetParametersAsVector()
{
    ListType params(this->GetNumberOfParameters(),0);

    params[0] = this->GetOrientationTheta();
    params[1] = this->GetOrientationPhi();
    params[2] = this->GetPerpendicularAngle();
    
    if (m_EstimateDiffusivities)
    {
        params[3] = this->GetAxialDiffusivity() - this->GetRadialDiffusivity1();
        params[4] = this->GetRadialDiffusivity1() - this->GetRadialDiffusivity2();
        params[5] = this->GetRadialDiffusivity2();
    }

    if (this->GetUseBoundedOptimization())
        this->UnboundParameters(params);

    return params;
}

TensorCompartment::ListType TensorCompartment::GetParameterLowerBounds()
{
    ListType lowerBounds(this->GetNumberOfParameters(),m_ZeroLowerBound);
    
    if (m_EstimateDiffusivities)
        lowerBounds[5] = m_DiffusivityLowerBound;
    
    return lowerBounds;
}

TensorCompartment::ListType TensorCompartment::GetParameterUpperBounds()
{
    ListType upperBounds(this->GetNumberOfParameters(),0);

    upperBounds[0] = m_PolarAngleUpperBound;
    upperBounds[1] = m_AzimuthAngleUpperBound;
    upperBounds[2] = m_AzimuthAngleUpperBound;

    if (m_EstimateDiffusivities)
    {
        upperBounds[3] = m_DiffusivityUpperBound;
        upperBounds[4] = m_DiffusivityUpperBound;
        upperBounds[5] = m_RadialDiffusivityUpperBound;
    }

    return upperBounds;
}
    
TensorCompartment::ListType TensorCompartment::BoundParameters(const ListType &params)
{
    ListType boundedParams(params);
    
    double inputSign = 1;
    boundedParams[0] = mcm_utilities::ComputeBoundedValue(params[0], inputSign, m_ZeroLowerBound, m_PolarAngleUpperBound);
    this->SetBoundedSignVectorValue(0,inputSign);
    boundedParams[1] = mcm_utilities::ComputeBoundedValue(params[1], inputSign, m_ZeroLowerBound, m_AzimuthAngleUpperBound);
    this->SetBoundedSignVectorValue(1,inputSign);
    boundedParams[2] = mcm_utilities::ComputeBoundedValue(params[2], inputSign, m_ZeroLowerBound, m_AzimuthAngleUpperBound);
    this->SetBoundedSignVectorValue(2,inputSign);

    if (m_EstimateDiffusivities)
    {
        boundedParams[3] = mcm_utilities::ComputeBoundedValue(params[3], inputSign, m_ZeroLowerBound, m_DiffusivityUpperBound);
        this->SetBoundedSignVectorValue(3,inputSign);
        boundedParams[4] = mcm_utilities::ComputeBoundedValue(params[4], inputSign, m_ZeroLowerBound, m_DiffusivityUpperBound);
        this->SetBoundedSignVectorValue(4,inputSign);
        boundedParams[5] = mcm_utilities::ComputeBoundedValue(params[5], inputSign, m_DiffusivityLowerBound, m_RadialDiffusivityUpperBound);
        this->SetBoundedSignVectorValue(5,inputSign);
    }

    return boundedParams;
}

void TensorCompartment::UnboundParameters(ListType &params)
{
    params[0] = mcm_utilities::UnboundValue(params[0], m_ZeroLowerBound, m_PolarAngleUpperBound);
    params[1] = mcm_utilities::UnboundValue(params[1], m_ZeroLowerBound, m_AzimuthAngleUpperBound);
    params[2] = mcm_utilities::UnboundValue(params[2], m_ZeroLowerBound, m_AzimuthAngleUpperBound);
    
    if (m_EstimateDiffusivities)
    {
        params[3] = mcm_utilities::UnboundValue(params[3], m_ZeroLowerBound, m_DiffusivityUpperBound);
        params[4] = mcm_utilities::UnboundValue(params[4], m_ZeroLowerBound, m_DiffusivityUpperBound);
        params[5] = mcm_utilities::UnboundValue(params[5], m_DiffusivityLowerBound, m_RadialDiffusivityUpperBound);
    }
}

void TensorCompartment::SetEstimateDiffusivities(bool arg)
{
    if (m_EstimateDiffusivities == arg)
        return;

    m_EstimateDiffusivities = arg;
    m_ChangedConstraints = true;
}

void TensorCompartment::SetCompartmentVector(ModelOutputVectorType &compartmentVector)
{
    if (compartmentVector.GetSize() != this->GetCompartmentSize())
        itkExceptionMacro("The input vector size does not match the size of the compartment");

    anima::GetTensorFromVectorRepresentation(compartmentVector,m_WorkVnlMatrix1,m_SpaceDimension);
    m_DiffusionTensor = m_WorkVnlMatrix1;

    m_TensorDeterminant = vnl_determinant <double> (m_WorkVnlMatrix1);

    if (m_TensorDeterminant > 0)
    {
        anima::GetTensorPower(m_WorkVnlMatrix1, m_WorkVnlMatrix2, -1.0);
        m_InverseDiffusionTensor = m_WorkVnlMatrix2;
    }
    else
    {
        m_TensorDeterminant = 0;
        m_InverseDiffusionTensor.Fill(0);
    }

    m_ModifiedTensor = false;
    m_ModifiedAngles = false;
    m_UpdatedCompartment = true;
}

void TensorCompartment::UpdateParametersFromTensor()
{
    if (!m_UpdatedCompartment)
        return;
    
    vnl_matrix <double> eVecs(m_SpaceDimension,m_SpaceDimension,0);
    vnl_diag_matrix <double> eVals(m_SpaceDimension);

    itk::SymmetricEigenAnalysis < Matrix3DType,vnl_diag_matrix <double>,vnl_matrix <double> > eigSys(m_SpaceDimension);

    eigSys.ComputeEigenValuesAndVectors(m_DiffusionTensor,eVals,eVecs);
    
    Vector3DType e1, sphDir;
    
    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        e1[i] = eVecs(2,i);

    anima::TransformCartesianToSphericalCoordinates(e1,sphDir);
    
    double theta = sphDir[0];
    double phi = sphDir[1];

    double alpha = 0;
    if ((theta != 0)&&(theta != M_PI))
        alpha = std::atan2(eVecs(1,2),- eVecs(0,2));
    else
    {
        // Use two first components of 2nd eigenvector
        // Infinite number of solutions for phi from main eigenvector, set it as alpha - pi/2
        alpha = std::atan2(1 - eVecs(1,0), eVecs(1,1));
        phi = alpha - M_PI / 2.0;
    }

    if (alpha < 0)
        alpha += 2.0 * M_PI;

    if (phi < 0)
        phi += 2.0 * M_PI;

    m_UpdatedCompartment = false;

    this->SetOrientationTheta(theta);
    this->SetOrientationPhi(phi);
    this->SetPerpendicularAngle(alpha);

    double axialDiff = std::max(m_DiffusivityLowerBound,std::min(eVals[2],m_DiffusivityUpperBound));
    double radialDiff1 = std::max(m_DiffusivityLowerBound,std::min(eVals[1],m_DiffusivityUpperBound));
    double radialDiff2 = std::max(m_DiffusivityLowerBound,std::min(eVals[0],m_DiffusivityUpperBound));

    this->SetAxialDiffusivity(axialDiff);
    this->SetRadialDiffusivity1(radialDiff1);
    this->SetRadialDiffusivity2(radialDiff2);
}

unsigned int TensorCompartment::GetCompartmentSize()
{
    return 6;
}

unsigned int TensorCompartment::GetNumberOfParameters()
{
    if (!m_ChangedConstraints)
        return m_NumberOfParameters;

    m_NumberOfParameters = this->GetCompartmentSize();

    if (!m_EstimateDiffusivities)
        m_NumberOfParameters -= 3;

    m_ChangedConstraints = false;
    return m_NumberOfParameters;
}

TensorCompartment::ModelOutputVectorType &TensorCompartment::GetCompartmentVector()
{
    if (m_CompartmentVector.GetSize() != this->GetCompartmentSize())
        m_CompartmentVector.SetSize(this->GetCompartmentSize());

    this->UpdateDiffusionTensor();

    anima::GetVectorRepresentation(m_DiffusionTensor.GetVnlMatrix().as_matrix(),m_CompartmentVector,this->GetCompartmentSize());
    return m_CompartmentVector;
}

void TensorCompartment::Reorient(vnl_matrix <double> &orientationMatrix, bool affineTransform)
{
    this->UpdateDiffusionTensor();

    m_WorkVnlMatrix1 = m_DiffusionTensor.GetVnlMatrix();

    if (!affineTransform)
        anima::RotateSymmetricMatrix(m_WorkVnlMatrix1,orientationMatrix,m_WorkVnlMatrix2);
    else
    {
        // PPD re-orientation
        vnl_matrix <double> ppdOrientationMatrix(3,3);
        typedef itk::Matrix <double, 3, 3> EigVecMatrixType;
        typedef vnl_vector_fixed <double,3> EigValVectorType;
        itk::SymmetricEigenAnalysis < EigVecMatrixType, EigValVectorType, EigVecMatrixType> eigenComputer(3);
        EigVecMatrixType eigVecs;
        EigValVectorType eigVals;

        eigenComputer.ComputeEigenValuesAndVectors(m_WorkVnlMatrix1,eigVals,eigVecs);
        anima::ExtractPPDRotationFromJacobianMatrix(orientationMatrix,ppdOrientationMatrix,eigVecs);
        anima::RotateSymmetricMatrix(m_WorkVnlMatrix1,ppdOrientationMatrix,m_WorkVnlMatrix2);
    }

    ModelOutputVectorType tmpVec(this->GetCompartmentSize());
    anima::GetVectorRepresentation(m_WorkVnlMatrix2,tmpVec);

    this->SetCompartmentVector(tmpVec);
}

const TensorCompartment::Matrix3DType &TensorCompartment::GetDiffusionTensor()
{
    this->UpdateDiffusionTensor();
    return m_DiffusionTensor;
}

void TensorCompartment::UpdateDiffusionTensor()
{
    if (!m_ModifiedTensor)
        return;

    this->UpdateAngleConfiguration();
    
    m_DiffusionTensor.Fill(0.0);
    
    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        m_DiffusionTensor(i,i) = this->GetRadialDiffusivity2();
    
    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        for (unsigned int j = i;j < m_SpaceDimension;++j)
        {
            m_DiffusionTensor(i,j) += m_EigenVector1[i] * m_EigenVector1[j] * (this->GetAxialDiffusivity() - this->GetRadialDiffusivity2());
            m_DiffusionTensor(i,j) += m_EigenVector2[i] * m_EigenVector2[j] * (this->GetRadialDiffusivity1() - this->GetRadialDiffusivity2());
            
            if (i != j)
                m_DiffusionTensor(j,i) = m_DiffusionTensor(i,j);
        }
    
    m_TensorDeterminant = this->GetAxialDiffusivity() * this->GetRadialDiffusivity1() * this->GetRadialDiffusivity2();
    
    m_WorkVnlMatrix1.set_size(m_SpaceDimension,m_SpaceDimension);
    m_WorkVnlMatrix1 = m_DiffusionTensor.GetVnlMatrix();
    anima::GetTensorPower(m_WorkVnlMatrix1, m_WorkVnlMatrix2, -1.0);
    m_InverseDiffusionTensor = m_WorkVnlMatrix2;
    
    m_ModifiedTensor = false;
}
    
void TensorCompartment::UpdateAngleConfiguration()
{
    if (!m_ModifiedAngles)
        return;

    m_CosTheta = std::cos(this->GetOrientationTheta());
    m_SinTheta = std::sin(this->GetOrientationTheta());
    m_CosPhi = std::cos(this->GetOrientationPhi());
    m_SinPhi = std::sin(this->GetOrientationPhi());
    m_CosAlpha = std::cos(this->GetPerpendicularAngle());
    m_SinAlpha = std::sin(this->GetPerpendicularAngle());
    
    m_EigenVector1[0] = m_SinTheta * m_CosPhi;
    m_EigenVector1[1] = m_SinTheta * m_SinPhi;
    m_EigenVector1[2] = m_CosTheta;
    
    m_EigenVector2[0] = -1.0 * (m_CosTheta * m_CosPhi * m_SinAlpha + m_SinPhi * m_CosAlpha);
    m_EigenVector2[1] = m_CosPhi * m_CosAlpha - m_CosTheta * m_SinPhi * m_SinAlpha;
    m_EigenVector2[2] = m_SinTheta * m_SinAlpha;

    m_ModifiedAngles = false;
}

} //end namespace anima
