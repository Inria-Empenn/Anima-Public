#include <animaBaseCompartment.h>

#include <animaVectorOperations.h>
#include <animaHyperbolicFunctions.h>

namespace anima
{

const double BaseCompartment::m_ZeroLowerBound = 0.0;
const double BaseCompartment::m_DiffusivityLowerBound = 1e-5;
const double BaseCompartment::m_AxialDiffusivityAddonLowerBound = 5.0e-4;
const double BaseCompartment::m_PolarAngleUpperBound = M_PI;
const double BaseCompartment::m_AzimuthAngleUpperBound = 2.0 * M_PI;
const double BaseCompartment::m_DiffusivityUpperBound = 3e-3;
const double BaseCompartment::m_RadialDiffusivityUpperBound = 1e-3;
const double BaseCompartment::m_DefaultConcentrationUpperBound = 20;
const double BaseCompartment::m_WatsonKappaUpperBound = 128.0;
const double BaseCompartment::m_Epsilon = 1.0e-2;
const double BaseCompartment::m_FractionUpperBound = 1.0;
const double BaseCompartment::m_TissuesRadiusLowerBound = 1.0e-7;
const double BaseCompartment::m_TissuesRadiusUpperBound = 40.1e-6;

double BaseCompartment::GetPredictedSignal(double smallDelta, double largeDelta, double gradientStrength, const Vector3DType &gradient)
{
    double ftDiffusionProfile = GetFourierTransformedDiffusionProfile(smallDelta,largeDelta, gradientStrength, gradient);

    return std::abs(ftDiffusionProfile);
}

bool BaseCompartment::IsEqual(Self *rhs, double tolerance)
{
    if (this->GetTensorCompatible() && rhs->GetTensorCompatible())
    {
        // Compare tensor representations, easier and probably faster
        Matrix3DType lhsTensor = this->GetDiffusionTensor();
        Matrix3DType rhsTensor = rhs->GetDiffusionTensor();

        for (unsigned int i = 0;i < Matrix3DType::RowDimensions;++i)
            for (unsigned int j = i;j < Matrix3DType::ColumnDimensions;++j)
            {
                if (std::abs(lhsTensor(i,j) - rhsTensor(i,j)) > tolerance)
                    return false;
            }

        return true;
    }

    if (std::abs(this->GetAxialDiffusivity() - rhs->GetAxialDiffusivity()) > tolerance)
        return false;

    Vector3DType orientationLHS(0.0);
    anima::TransformSphericalToCartesianCoordinates(this->GetOrientationTheta(),this->GetOrientationPhi(),1.0,orientationLHS);
    Vector3DType orientationRHS(0.0);
    anima::TransformSphericalToCartesianCoordinates(rhs->GetOrientationTheta(),rhs->GetOrientationPhi(),1.0,orientationRHS);

    if (std::abs(std::abs(anima::ComputeScalarProduct(orientationLHS,orientationRHS)) - 1.0) > tolerance)
        return false;
    if (std::abs(this->GetPerpendicularAngle() - rhs->GetPerpendicularAngle()) > tolerance)
        return false;
    if (std::abs(this->GetRadialDiffusivity1() - rhs->GetRadialDiffusivity1()) > tolerance)
        return false;
    if (std::abs(this->GetRadialDiffusivity2() - rhs->GetRadialDiffusivity2()) > tolerance)
        return false;
    if (std::abs(this->GetOrientationConcentration() - rhs->GetOrientationConcentration()) > tolerance)
        return false;
    if (std::abs(this->GetExtraAxonalFraction() - rhs->GetExtraAxonalFraction()) > tolerance)
        return false;
    if (std::abs(this->GetTissueRadius() - rhs->GetTissueRadius()) > tolerance)
        return false;

    return true;
}

void BaseCompartment::CopyFromOther(Self *rhs)
{
    this->SetOrientationTheta(rhs->GetOrientationTheta());
    this->SetOrientationPhi(rhs->GetOrientationPhi());
    this->SetPerpendicularAngle(rhs->GetPerpendicularAngle());
    this->SetAxialDiffusivity(rhs->GetAxialDiffusivity());
    this->SetRadialDiffusivity1(rhs->GetRadialDiffusivity1());
    this->SetRadialDiffusivity2(rhs->GetRadialDiffusivity2());
    this->SetOrientationConcentration(rhs->GetOrientationConcentration());
    this->SetExtraAxonalFraction(rhs->GetExtraAxonalFraction());
    this->SetTissueRadius(rhs->GetTissueRadius());
}

void BaseCompartment::Reorient(vnl_matrix <double> &orientationMatrix, bool affineTransform)
{
    Vector3DType tmpDir;
    Vector3DType tmpRotatedDir;

    anima::TransformSphericalToCartesianCoordinates(this->GetOrientationTheta(),this->GetOrientationPhi(),1.0,tmpDir);

    unsigned int imageDimension = tmpDir.size();
    for (unsigned int k = 0;k < imageDimension;++k)
    {
        tmpRotatedDir[k] = 0;

        for (unsigned int l = 0;l < imageDimension;++l)
            tmpRotatedDir[k] += orientationMatrix(l,k) * tmpDir[l];
    }

    anima::TransformCartesianToSphericalCoordinates(tmpRotatedDir,tmpRotatedDir);
    this->SetOrientationTheta(tmpRotatedDir[0]);
    this->SetOrientationPhi(tmpRotatedDir[1]);
}

const BaseCompartment::Matrix3DType &BaseCompartment::GetDiffusionTensor()
{
    throw itk::ExceptionObject(__FILE__,__LINE__,"This compartment type does not support diffusion tensor export",ITK_LOCATION);
}

double BaseCompartment::GetFractionalAnisotropy()
{
    throw itk::ExceptionObject(__FILE__,__LINE__,"This compartment type does not support fractional anisotropy computation",ITK_LOCATION);
}

} // end namespace anima
