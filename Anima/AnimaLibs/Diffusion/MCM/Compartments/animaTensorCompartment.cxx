#include <animaTensorCompartment.h>

#include <itkSymmetricEigenAnalysis.h>
#include <animaMCMConstants.h>

namespace anima
{

double TensorCompartment::GetFourierTransformedDiffusionProfile(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient)
{
    this->UpdateDiffusionTensor();

    double quadForm = 0;

    for (unsigned int i = 0;i < m_SpaceDimension;++i)
    {
        quadForm += m_DiffusionTensor(i,i) * gradient[i] * gradient[i];
        for (unsigned int j = i+1;j < m_SpaceDimension;++j)
            quadForm += 2 * m_DiffusionTensor(i,j) * gradient[i] * gradient[j];
    }

    double bValue = anima::GetBValueFromAcquisitionParameters(smallDelta, bigDelta, gradientStrength);

    return std::exp(- bValue * quadForm);
}

TensorCompartment::ListType &TensorCompartment::GetSignalAttenuationJacobian(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient)
{
    this->UpdateDiffusionTensor();

    m_JacobianVector.resize(this->GetNumberOfParameters());
    
    double signalAttenuation = this->GetFourierTransformedDiffusionProfile(smallDelta, bigDelta, gradientStrength, gradient);
    double innerProd1 = anima::ComputeScalarProduct(gradient, m_EigenVector1);
    double innerProd2 = anima::ComputeScalarProduct(gradient, m_EigenVector2);
    
    double DgTe1DTheta = m_CosTheta * (gradient[0] * m_CosPhi + gradient[1] * m_SinPhi) - gradient[2] * m_SinTheta;
    double DgTe1DPhi = m_SinTheta * (gradient[1] * m_CosPhi - gradient[0] * m_SinPhi);
    
    double DgTe2DTheta = m_SinAlpha * innerProd1;
    double DgTe2DPhi = gradient[0] * (m_CosTheta * m_SinPhi * m_SinAlpha - m_CosPhi * m_CosAlpha) - gradient[1] * (m_SinPhi * m_CosAlpha + m_CosTheta * m_CosPhi * m_SinAlpha);
    double DgTe2DAlpha = gradient[0] * (m_SinPhi * m_SinAlpha - m_CosTheta * m_CosPhi * m_CosAlpha) - gradient[1] * (m_CosPhi * m_SinAlpha + m_CosTheta * m_SinPhi * m_CosAlpha) + gradient[2] * m_SinTheta * m_CosAlpha;
    
    double diffAxialRadial2 = this->GetAxialDiffusivity() - this->GetRadialDiffusivity2();
    double diffRadialDiffusivities = this->GetRadialDiffusivity1() - this->GetRadialDiffusivity2();

    // Derivative w.r.t. theta
    double bValue = anima::GetBValueFromAcquisitionParameters(smallDelta, bigDelta, gradientStrength);
    m_JacobianVector[0] = -2.0 * bValue * (diffAxialRadial2 * innerProd1 * DgTe1DTheta + diffRadialDiffusivities * innerProd2 * DgTe2DTheta) * signalAttenuation;
    
    // Derivative w.r.t. phi
    m_JacobianVector[1] = -2.0 * bValue * (diffAxialRadial2 * innerProd1 * DgTe1DPhi + diffRadialDiffusivities * innerProd2 * DgTe2DPhi) * signalAttenuation;
    
    // Derivative w.r.t. alpha
    m_JacobianVector[2] = -2.0 * bValue * diffRadialDiffusivities * innerProd2 * DgTe2DAlpha * signalAttenuation;

    if (m_EstimateDiffusivities)
    {
        // Derivative w.r.t. to d1
        m_JacobianVector[3] = - bValue * innerProd1 * innerProd1 * signalAttenuation;
        
        // Derivative w.r.t. to d2
        m_JacobianVector[4] = - bValue * (innerProd1 * innerProd1 + innerProd2 * innerProd2) * signalAttenuation;
        
        // Derivative w.r.t. to d3
        m_JacobianVector[5] = - bValue * signalAttenuation;
    }
    
    return m_JacobianVector;
}

double TensorCompartment::GetLogDiffusionProfile(const Vector3DType &sample)
{
    this->UpdateInverseDiffusionTensor();

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

    this->SetOrientationTheta(params[0]);
    this->SetOrientationPhi(params[1]);
    this->SetPerpendicularAngle(params[2]);

    if (m_EstimateDiffusivities)
    {
        this->SetAxialDiffusivity(params[3] + params[4] + params[5]);
        this->SetRadialDiffusivity1(params[4] + params[5]);
        this->SetRadialDiffusivity2(params[5]);
    }
}

TensorCompartment::ListType &TensorCompartment::GetParametersAsVector()
{
    m_ParametersVector.resize(this->GetNumberOfParameters());

    m_ParametersVector[0] = this->GetOrientationTheta();
    m_ParametersVector[1] = this->GetOrientationPhi();
    m_ParametersVector[2] = this->GetPerpendicularAngle();
    
    if (m_EstimateDiffusivities)
    {
        m_ParametersVector[3] = this->GetAxialDiffusivity() - this->GetRadialDiffusivity1();
        m_ParametersVector[4] = this->GetRadialDiffusivity1() - this->GetRadialDiffusivity2();
        m_ParametersVector[5] = this->GetRadialDiffusivity2();
    }

    return m_ParametersVector;
}

TensorCompartment::ListType &TensorCompartment::GetParameterLowerBounds()
{
    m_ParametersLowerBoundsVector.resize(this->GetNumberOfParameters());
    std::fill(m_ParametersLowerBoundsVector.begin(),m_ParametersLowerBoundsVector.end(),anima::MCMZeroLowerBound);
    
    if (m_EstimateDiffusivities)
    {
        m_ParametersLowerBoundsVector[3] = anima::MCMAxialDiffusivityAddonLowerBound;
        m_ParametersLowerBoundsVector[5] = anima::MCMDiffusivityLowerBound;
    }
    
    return m_ParametersLowerBoundsVector;
}

TensorCompartment::ListType &TensorCompartment::GetParameterUpperBounds()
{
    m_ParametersUpperBoundsVector.resize(this->GetNumberOfParameters());

    m_ParametersUpperBoundsVector[0] = anima::MCMPolarAngleUpperBound;
    m_ParametersUpperBoundsVector[1] = anima::MCMAzimuthAngleUpperBound;
    m_ParametersUpperBoundsVector[2] = anima::MCMAzimuthAngleUpperBound;

    if (m_EstimateDiffusivities)
    {
        m_ParametersUpperBoundsVector[3] = anima::MCMDiffusivityUpperBound;
        m_ParametersUpperBoundsVector[4] = anima::MCMDiffusivityUpperBound;
        m_ParametersUpperBoundsVector[5] = anima::MCMRadialDiffusivityUpperBound;
    }

    return m_ParametersUpperBoundsVector;
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

    if (m_TensorDeterminant <= 0)
        m_TensorDeterminant = 0;

    m_ModifiedTensor = false;
    m_ModifiedAngles = false;
    m_UpdateInverseTensor = true;
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

    double axialDiff = std::max(anima::MCMDiffusivityLowerBound,std::min(eVals[2],anima::MCMDiffusivityUpperBound));
    double radialDiff1 = std::max(anima::MCMDiffusivityLowerBound,std::min(eVals[1],anima::MCMDiffusivityUpperBound));
    double radialDiff2 = std::max(anima::MCMDiffusivityLowerBound,std::min(eVals[0],anima::MCMDiffusivityUpperBound));

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

    m_WorkVnlMatrix1 = m_DiffusionTensor.GetVnlMatrix().as_matrix();

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
    double axialDiff = this->GetAxialDiffusivity();
    double radialDiff1 = this->GetRadialDiffusivity1();
    double radialDiff2 = this->GetRadialDiffusivity2();
    
    m_DiffusionTensor.Fill(0.0);
    
    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        m_DiffusionTensor(i,i) = radialDiff2;
    
    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        for (unsigned int j = i;j < m_SpaceDimension;++j)
        {
            m_DiffusionTensor(i,j) += m_EigenVector1[i] * m_EigenVector1[j] * (axialDiff - radialDiff2);
            m_DiffusionTensor(i,j) += m_EigenVector2[i] * m_EigenVector2[j] * (radialDiff1 - radialDiff2);
            
            if (i != j)
                m_DiffusionTensor(j,i) = m_DiffusionTensor(i,j);
        }
    
    m_TensorDeterminant = axialDiff * radialDiff1 * radialDiff2;
    
    m_UpdateInverseTensor = true;
    m_ModifiedTensor = false;
}

void TensorCompartment::UpdateInverseDiffusionTensor()
{
    this->UpdateDiffusionTensor();

    if (!m_UpdateInverseTensor)
        return;

    if (m_TensorDeterminant == 0.0)
    {
        m_InverseDiffusionTensor.Fill(0.0);
        m_UpdateInverseTensor = false;
        return;
    }

    m_WorkVnlMatrix1.set_size(m_SpaceDimension,m_SpaceDimension);
    m_WorkVnlMatrix1 = m_DiffusionTensor.GetVnlMatrix().as_matrix();
    m_leCalculator->GetTensorPower(m_WorkVnlMatrix1, m_WorkVnlMatrix2, -1.0);
    m_InverseDiffusionTensor = m_WorkVnlMatrix2;

    m_UpdateInverseTensor = false;
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

double TensorCompartment::GetApparentFractionalAnisotropy()
{
    double l1 = this->GetAxialDiffusivity();
    double l2 = this->GetRadialDiffusivity1();
    double l3 = this->GetRadialDiffusivity2();

    double numFA = std::sqrt ((l1 - l2) * (l1 - l2) + (l2 - l3) * (l2 - l3) + (l3 - l1) * (l3 - l1));
    double denomFA = std::sqrt (l1 * l1 + l2 * l2 + l3 * l3);

    double fa = 0;
    if (denomFA != 0)
        fa = std::sqrt(0.5) * (numFA / denomFA);

    return fa;
}

double TensorCompartment::GetApparentMeanDiffusivity()
{
    double l1 = this->GetAxialDiffusivity();
    double l2 = this->GetRadialDiffusivity1();
    double l3 = this->GetRadialDiffusivity2();

    return (l1 + l2 + l3) / 3.0;
}

double TensorCompartment::GetApparentParallelDiffusivity()
{
    return this->GetAxialDiffusivity();
}

double TensorCompartment::GetApparentPerpendicularDiffusivity()
{
    return (this->GetRadialDiffusivity1() + this->GetRadialDiffusivity2()) / 2.0;
}

} //end namespace anima
