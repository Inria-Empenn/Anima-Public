#include <animaBaseCompartment.h>

#include <animaVectorOperations.h>
#include <animaHyperbolicFunctions.h>

namespace anima
{

double BaseCompartment::GetPredictedSignal(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient)
{
    double ftDiffusionProfile = GetFourierTransformedDiffusionProfile(smallDelta,bigDelta, gradientStrength, gradient);

    return std::abs(ftDiffusionProfile);
}

bool BaseCompartment::IsEqual(Self *rhs, double tolerance, double absoluteTolerance)
{
    if (this->GetTensorCompatible() && rhs->GetTensorCompatible())
    {
        // Compare tensor representations, easier and probably faster
        Matrix3DType lhsTensor = this->GetDiffusionTensor();
        Matrix3DType rhsTensor = rhs->GetDiffusionTensor();

        for (unsigned int i = 0;i < Matrix3DType::RowDimensions;++i)
        {
            for (unsigned int j = i;j < Matrix3DType::ColumnDimensions;++j)
            {
                double diffValue = std::abs(lhsTensor(i,j) - rhsTensor(i,j));
                double denom = std::max(std::abs(lhsTensor(i,j)),std::abs(rhsTensor(i,j)));
                if (denom == 0)
                    continue;
                if ((diffValue / denom > tolerance) && (diffValue > absoluteTolerance))
                    return false;
            }
        }

        return true;
    }

    double denomValue = std::max(this->GetAxialDiffusivity(),rhs->GetAxialDiffusivity());
    double diffValue = std::abs(this->GetAxialDiffusivity() - rhs->GetAxialDiffusivity());
    if ((diffValue / denomValue > tolerance) && (diffValue > absoluteTolerance))
        return false;

    Vector3DType orientationLHS(0.0);
    anima::TransformSphericalToCartesianCoordinates(this->GetOrientationTheta(),this->GetOrientationPhi(),1.0,orientationLHS);
    Vector3DType orientationRHS(0.0);
    anima::TransformSphericalToCartesianCoordinates(rhs->GetOrientationTheta(),rhs->GetOrientationPhi(),1.0,orientationRHS);

    if (std::abs(std::abs(anima::ComputeScalarProduct(orientationLHS,orientationRHS)) - 1.0) > tolerance)
        return false;

    denomValue = std::max(this->GetPerpendicularAngle(),rhs->GetPerpendicularAngle());
    diffValue = std::abs(this->GetPerpendicularAngle() - rhs->GetPerpendicularAngle());
    if ((diffValue / denomValue > tolerance) && (diffValue > absoluteTolerance))
        return false;

    denomValue = std::max(this->GetRadialDiffusivity1(),rhs->GetRadialDiffusivity1());
    diffValue = std::abs(this->GetRadialDiffusivity1() - rhs->GetRadialDiffusivity1());
    if ((diffValue / denomValue > tolerance) && (diffValue > absoluteTolerance))
        return false;

    denomValue = std::max(this->GetRadialDiffusivity2(),rhs->GetRadialDiffusivity2());
    diffValue = std::abs(this->GetRadialDiffusivity2() - rhs->GetRadialDiffusivity2());
    if ((diffValue / denomValue > tolerance) && (diffValue > absoluteTolerance))
        return false;

    denomValue = std::max(this->GetOrientationConcentration(),rhs->GetOrientationConcentration());
    diffValue = std::abs(this->GetOrientationConcentration() - rhs->GetOrientationConcentration());
    if ((diffValue / denomValue > tolerance) && (diffValue > absoluteTolerance))
        return false;

    denomValue = std::max(this->GetExtraAxonalFraction(),rhs->GetExtraAxonalFraction());
    diffValue = std::abs(this->GetExtraAxonalFraction() - rhs->GetExtraAxonalFraction());
    if ((diffValue / denomValue > tolerance) && (diffValue > absoluteTolerance))
        return false;

    denomValue = std::max(this->GetTissueRadius(),rhs->GetTissueRadius());
    diffValue = std::abs(this->GetTissueRadius() - rhs->GetTissueRadius());
    if ((diffValue / denomValue > tolerance) && (diffValue > absoluteTolerance))
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
