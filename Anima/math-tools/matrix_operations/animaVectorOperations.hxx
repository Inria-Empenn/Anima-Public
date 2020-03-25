#include "animaVectorOperations.h"
#include <cmath>

#include <animaLogarithmFunctions.h>

#include <vnl/vnl_diag_matrix.h>
#include <itkMacro.h>

namespace anima
{

/******* Main function GetMedian *******/
// Main
template <class VectorType> double GetMedian(const VectorType &data, const unsigned int NDimension)
{
    std::vector <double> array(NDimension, 0);
    
    for (unsigned int i = 0;i < NDimension;++i)
        array[i] = data[i];
    
    std::nth_element(array.begin(), array.begin() + NDimension / 2, array.end());
    
    double median = array[NDimension / 2];
    
    if (NDimension % 2 == 0)
    {
        median += *std::max_element(array.begin(), array.begin() + NDimension / 2);
        median /= 2.0;
    }
    
    return median;
}
// For itkVector
template <class ScalarType, unsigned int NDimension> double GetMedian(const itk::Vector <ScalarType,NDimension> &data)
{
    return GetMedian(data, NDimension);
}
// For itkVariableLengthVector
template <class ScalarType> double GetMedian(const itk::VariableLengthVector <ScalarType> &data)
{
    unsigned int NDimension = data.GetSize();
    return GetMedian(data, NDimension);
}
// For itkPoint
template <class ScalarType, unsigned int NDimension> double GetMedian(const itk::Point<ScalarType,NDimension> &data)
{
    return GetMedian(data, NDimension);
}
// For vnl_vector
template <class ScalarType> double GetMedian(const vnl_vector <ScalarType> &data)
{
    unsigned int NDimension = data.size();
    return GetMedian(data, NDimension);
}
// For vnl_vector_fixed
template <class ScalarType, unsigned int NDimension>
double GetMedian(const vnl_vector_fixed <ScalarType,NDimension> &data)
{
    return GetMedian(data, NDimension);
}
// For std::vector
template <class ScalarType> double GetMedian(const std::vector <ScalarType> &data)
{
    unsigned int NDimension = data.size();
    return GetMedian(data, NDimension);
}
/******************************************************/




/******* Main function ComputeEuclideanDistance *******/
// Main
template <class VectorType> double ComputeEuclideanDistance(const VectorType &x1, const VectorType &x2, const unsigned int NDimension)
{
    double resVal = 0;

    for (unsigned int i = 0;i < NDimension;++i)
    {
        resVal += (x1[i] - x2[i]) * (x1[i] - x2[i]);
    }

    return std::sqrt(resVal);
}
// For itkVector
template <class ScalarType, unsigned int NDimension> double ComputeEuclideanDistance(const itk::Vector <ScalarType,NDimension> &x1, const itk::Vector <ScalarType,NDimension> &x2)
{
    return ComputeEuclideanDistance(x1, x2, NDimension);
}
// For itkVariableLengthVector
template <class ScalarType> double ComputeEuclideanDistance(const itk::VariableLengthVector <ScalarType> &x1, const itk::VariableLengthVector <ScalarType> &x2)
{
    unsigned int NDimension = x1.GetSize();

    if (x2.GetSize() != NDimension)
        throw itk::ExceptionObject(__FILE__, __LINE__,"Euclidean distance can only be computed on same size vectors.",ITK_LOCATION);

    return ComputeEuclideanDistance(x1, x2, NDimension);
}
// For itkPoint
template <class ScalarType, unsigned int NDimension> double ComputeEuclideanDistance(const itk::Point<ScalarType,NDimension> &x1, const itk::Point<ScalarType,NDimension> &x2)
{
    return ComputeEuclideanDistance(x1, x2, NDimension);
}
// For vnl_vector
template <class ScalarType> double ComputeEuclideanDistance(const vnl_vector <ScalarType> &x1, const vnl_vector <ScalarType> &x2)
{
    unsigned int NDimension = x1.size();

    if (x2.size() != NDimension)
        throw itk::ExceptionObject(__FILE__, __LINE__,"Euclidean distance can only be computed on same size vectors.",ITK_LOCATION);

    return ComputeEuclideanDistance(x1, x2, NDimension);
}
// For vnl_vector_fixed
template <class ScalarType, unsigned int NDimension>
double ComputeEuclideanDistance(const vnl_vector_fixed <ScalarType,NDimension> &x1, const vnl_vector_fixed <ScalarType,NDimension> &x2)
{
    return ComputeEuclideanDistance(x1,x2,NDimension);
}
// For std::vector
template <class ScalarType> double ComputeEuclideanDistance(const std::vector <ScalarType> &x1, const std::vector <ScalarType> &x2)
{
    unsigned int NDimension = x1.size();

    if (x2.size() != NDimension)
        throw itk::ExceptionObject(__FILE__, __LINE__,"Euclidean distance can only be computed on same size vectors.",ITK_LOCATION);

    return ComputeEuclideanDistance(x1, x2, NDimension);
}
/******************************************************/

// Function ComputePointToSetDistance
template <class VectorType> double ComputePointToSetDistance(const VectorType &x, const std::vector <VectorType> &s)
{
    unsigned int NDimension = s.size();

    double resVal = ComputeEuclideanDistance(x, s[0]);

    for (unsigned int i = 1;i < NDimension;++i)
    {
        double tmpVal = ComputeEuclideanDistance(x, s[i]);

        if (tmpVal < resVal)
            resVal = tmpVal;
    }

    return resVal;
}

// Function ComputeDirectedHausdorffDistance
template <class VectorType> double ComputeDirectedHausdorffDistance(const std::vector <VectorType> &s1, const std::vector <VectorType> &s2)
{
    unsigned int NDimension = s1.size();

    double resVal = ComputePointToSetDistance(s1[0], s2);

    for (unsigned int i = 1;i < NDimension;++i)
    {
        double tmpVal = ComputePointToSetDistance(s1[i], s2);

        if (tmpVal > resVal)
            resVal = tmpVal;
    }

    return resVal;
}

// Function ComputeHausdorffDistance
template <class VectorType> double ComputeHausdorffDistance(const std::vector <VectorType> &s1, const std::vector <VectorType> &s2)
{
    double  forward = ComputeDirectedHausdorffDistance(s1, s2);
    double backward = ComputeDirectedHausdorffDistance(s2, s1);

    return std::max(forward, backward);
}

// Function ComputeModifiedDirectedHausdorffDistance
template <class VectorType> double ComputeModifiedDirectedHausdorffDistance(const std::vector <VectorType> &s1, const std::vector <VectorType> &s2)
{
    unsigned int NDimension = s1.size();

    double resVal = ComputePointToSetDistance(s1[0], s2);

    for (unsigned int i = 1;i < NDimension;++i)
        resVal += ComputePointToSetDistance(s1[i], s2);

    return resVal / NDimension;
}

// Function ComputeModifiedHausdorffDistance
template <class VectorType> double ComputeModifiedHausdorffDistance(const std::vector <VectorType> &s1, const std::vector <VectorType> &s2)
{
    double  forward = ComputeModifiedDirectedHausdorffDistance(s1, s2);
    double backward = ComputeModifiedDirectedHausdorffDistance(s2, s1);

    return std::max(forward, backward);
}

/******* Function ExponentialSum *******/
// Main
template <class VectorType> double ExponentialSum(const VectorType &x, const unsigned int NDimension)
{
    double maxVal = x[0];

    for (unsigned int i = 1;i < NDimension;++i)
        if (x[i] > maxVal)
            maxVal = x[i];

    double sum = 0;
    for (unsigned int i = 0;i < NDimension;++i)
        sum += std::exp(x[i] - maxVal);

    if (sum < 1.0e-6)
        throw itk::ExceptionObject(__FILE__, __LINE__,"Exponential sum cannot be computed because requires to take the log of a non-positive value.",ITK_LOCATION);

    return maxVal + std::log(sum);
}
// For itkVector
template <class ScalarType, unsigned int NDimension> double ExponentialSum(const itk::Vector <ScalarType,NDimension> &x)
{
    return ExponentialSum(x, NDimension);
}
// For itkVariableLengthVector
template <class ScalarType> double ExponentialSum(const itk::VariableLengthVector <ScalarType> &x)
{
    unsigned int NDimension = x.GetSize();
    return ExponentialSum(x, NDimension);
}
// For itkPoint
template <class ScalarType, unsigned int NDimension> double ExponentialSum(const itk::Point <ScalarType,NDimension> &x)
{
    return ExponentialSum(x, NDimension);
}
// For vnl_vector
template <class ScalarType> double ExponentialSum(const vnl_vector <ScalarType> &x)
{
    unsigned int NDimension = x.size();
    return ExponentialSum(x, NDimension);
}
// For vnl_vector_fixed
template <class ScalarType, unsigned int NDimension> double ExponentialSum(const vnl_vector_fixed <ScalarType,NDimension> &x)
{
    return ExponentialSum(x, NDimension);
}
// For std::vector
template <class ScalarType> double ExponentialSum(const std::vector <ScalarType> &x)
{
    unsigned int NDimension = x.size();
    return ExponentialSum(x, NDimension);
}
/***************************************/

/******* Main function ComputeScalarProduct *******/
// Main
template <class VectorType> double ComputeScalarProduct(const VectorType &v1, const VectorType &v2, const unsigned int NDimension)
{
    double resVal = 0;

    for (unsigned int i = 0;i < NDimension;++i)
        resVal += v1[i] * v2[i];

    return resVal;
}
// For itkVector
template <class ScalarType, unsigned int NDimension> double ComputeScalarProduct(const itk::Vector <ScalarType,NDimension> &v1, const itk::Vector <ScalarType,NDimension> &v2)
{
    return ComputeScalarProduct(v1, v2, NDimension);
}
// For itkVariableLengthVector
template <class ScalarType> double ComputeScalarProduct(const itk::VariableLengthVector <ScalarType> &v1, const itk::VariableLengthVector <ScalarType> &v2)
{
    unsigned int NDimension = v1.GetSize();

    if (v2.GetSize() != NDimension)
        throw itk::ExceptionObject(__FILE__, __LINE__,"Scalar product can only be computed on same size vectors.",ITK_LOCATION);

    return ComputeScalarProduct(v1, v2, NDimension);
}
// For itkPoint
template <class ScalarType, unsigned int NDimension> double ComputeScalarProduct(const itk::Point <ScalarType,NDimension> &v1, const itk::Point <ScalarType,NDimension> &v2)
{
    return ComputeScalarProduct(v1, v2, NDimension);
}
// For vnl_vector
template <class ScalarType> double ComputeScalarProduct(const vnl_vector <ScalarType> &v1, const vnl_vector <ScalarType> &v2)
{
    unsigned int NDimension = v1.size();

    if (v2.size() != NDimension)
        throw itk::ExceptionObject(__FILE__, __LINE__,"Scalar product can only be computed on same size vectors.",ITK_LOCATION);

    return ComputeScalarProduct(v1, v2, NDimension);
}
// For vnl_vector_fixed
template <class ScalarType, unsigned int NDimension> double ComputeScalarProduct(const vnl_vector_fixed <ScalarType,NDimension> &v1, const vnl_vector_fixed <ScalarType,NDimension> &v2)
{
    return ComputeScalarProduct(v1, v2, NDimension);
}
// For std::vector
template <class ScalarType> double ComputeScalarProduct(const std::vector <ScalarType> &v1, const std::vector <ScalarType> &v2)
{
    unsigned int NDimension = v1.size();

    if (v2.size() != NDimension)
        throw itk::ExceptionObject(__FILE__, __LINE__,"Scalar product can only be computed on same size vectors.",ITK_LOCATION);

    return ComputeScalarProduct(v1, v2, NDimension);
}
/**************************************************/

/******* Function ComputeCrossProduct *******/
// Main
template <class VectorType> void ComputeCrossProduct(const VectorType &v1, const VectorType &v2, const unsigned int NDimension, VectorType &resVec)
{
    if (NDimension != 3)
        throw itk::ExceptionObject(__FILE__, __LINE__,"Cross product requires 3D vectors.",ITK_LOCATION);

    resVec[0] = v1[1] * v2[2] - v2[1] * v1[2];
    resVec[1] = v1[2] * v2[0] - v2[2] * v1[0];
    resVec[2] = v1[0] * v2[1] - v2[0] * v1[1];
}
// For itkVector
template <class ScalarType, unsigned int NDimension> void ComputeCrossProduct(const itk::Vector <ScalarType,NDimension> &v1, const itk::Vector <ScalarType,NDimension> &v2, itk::Vector <ScalarType,NDimension> &resVec)
{
    return ComputeCrossProduct(v1, v2, NDimension, resVec);
}
// For itkVariableLengthVector
template <class ScalarType> void ComputeCrossProduct(const itk::VariableLengthVector <ScalarType> &v1, const itk::VariableLengthVector <ScalarType> &v2, itk::VariableLengthVector <ScalarType> &resVec)
{
    unsigned int NDimension = v1.GetSize();

    if (v2.GetSize() != NDimension)
        throw itk::ExceptionObject(__FILE__, __LINE__,"Cross product can only be computed on same size vectors.",ITK_LOCATION);

    resVec.SetSize(NDimension);
    return ComputeCrossProduct(v1, v2, NDimension, resVec);
}
// For itkPoint
template <class ScalarType, unsigned int NDimension> void ComputeCrossProduct(const itk::Point <ScalarType,NDimension> &v1, const itk::Point <ScalarType,NDimension> &v2, itk::Point <ScalarType,NDimension> &resVec)
{
    return ComputeCrossProduct(v1, v2, NDimension, resVec);
}
// For vnl_vector
template <class ScalarType> void ComputeCrossProduct(const vnl_vector <ScalarType> &v1, const vnl_vector <ScalarType> &v2, vnl_vector <ScalarType> &resVec)
{
    unsigned int NDimension = v1.size();

    if (v2.size() != NDimension)
        throw itk::ExceptionObject(__FILE__, __LINE__,"Cross product can only be computed on same size vectors.",ITK_LOCATION);

    resVec.set_size(NDimension);
    return ComputeCrossProduct(v1, v2, NDimension, resVec);
}
// For vnl_vector_fixed
template <class ScalarType, unsigned int NDimension> void ComputeCrossProduct(const vnl_vector_fixed <ScalarType,NDimension> &v1, const vnl_vector_fixed <ScalarType,NDimension> &v2, vnl_vector_fixed <ScalarType,NDimension> &resVec)
{
    return ComputeCrossProduct(v1, v2, NDimension, resVec);
}
// For std::vector
template <class ScalarType> void ComputeCrossProduct(const std::vector <ScalarType> &v1, const std::vector <ScalarType> &v2, std::vector <ScalarType> &resVec)
{
    unsigned int NDimension = v1.size();

    if (v2.size() != NDimension)
        throw itk::ExceptionObject(__FILE__, __LINE__,"Cross product can only be computed on same size vectors.",ITK_LOCATION);

    resVec.resize(NDimension);
    return ComputeCrossProduct(v1, v2, NDimension, resVec);
}
/********************************************/

// Function ComputeNorm
template <class VectorType> double ComputeNorm(const VectorType &v)
{
    double n = ComputeScalarProduct(v, v);

    if (n < 1.0e-6)
        return 0;

    return std::sqrt(n);
}

/******* Function Normalize *******/
// Main
template <class VectorType> void Normalize(const VectorType &v, const unsigned int NDimension, VectorType &resVec)
{
    resVec = v;

    double n = ComputeNorm(v);

    if (n < 1.0e-6)
        return;

    for (unsigned int i = 0;i < NDimension;++i)
        resVec[i] /= n;
}
// For itkVector
template <class ScalarType, unsigned int NDimension> void Normalize(const itk::Vector <ScalarType,NDimension> &v, itk::Vector <ScalarType,NDimension> &resVec)
{
    Normalize(v, NDimension, resVec);
}
// For itkVariableLengthVector
template <class ScalarType> void Normalize(const itk::VariableLengthVector <ScalarType> &v, itk::VariableLengthVector <ScalarType> &resVec)
{
    unsigned int NDimension = v.GetSize();
    Normalize(v, NDimension, resVec);
}
// For itkPoint
template <class ScalarType, unsigned int NDimension> void Normalize(const itk::Point <ScalarType,NDimension> &v, itk::Point <ScalarType,NDimension> &resVec)
{
    Normalize(v, NDimension, resVec);
}
// For vnl_vector
template <class ScalarType> void Normalize(const vnl_vector <ScalarType> &v, vnl_vector <ScalarType> &resVec)
{
    unsigned int NDimension = v.size();
    Normalize(v, NDimension, resVec);
}
// For vnl_vector_fixed
template <class ScalarType, unsigned int NDimension> void Normalize(const vnl_vector_fixed <ScalarType, NDimension> &v, vnl_vector_fixed <ScalarType, NDimension> &resVec)
{
    Normalize(v, NDimension, resVec);
}
// For std::vector
template <class ScalarType> void Normalize(const std::vector <ScalarType> &v, std::vector <ScalarType> &resVec)
{
    unsigned int NDimension = v.size();
    Normalize(v, NDimension, resVec);
}
/**********************************/

/******* Householder vector functions *******/

// Main
template <class VectorType> void ComputeHouseholderVector(const VectorType &inputVector, VectorType &outputVector, double &beta, unsigned int dimension)
{
    double sigma = 0.0;
    outputVector[0] = 1.0;
    for (unsigned int i = 1;i < dimension;++i)
    {
        sigma += inputVector[i] * inputVector[i];
        outputVector[i] = inputVector[i];
    }

    if (sigma == 0.0)
        beta = 0.0;
    else
    {
        double mu = std::sqrt(inputVector[0] * inputVector[0] + sigma);
        if (inputVector[0] <= 0.0)
            outputVector[0] = inputVector[0] - mu;
        else
            outputVector[0] = - sigma / (inputVector[0] + mu);

        beta = 2 * outputVector[0] * outputVector[0] / (sigma + outputVector[0] * outputVector[0]);

        for (unsigned int i = 1;i < dimension;++i)
            outputVector[i] /= outputVector[0];

        outputVector[0] = 1.0;
    }
}

// for std vector
template <class ScalarType> void ComputeHouseholderVector(const std::vector <ScalarType> &inputVector, std::vector <ScalarType> &outputVector, double &beta)
{
    unsigned int dimension = inputVector.size();
    outputVector.resize(dimension);

    ComputeHouseholderVector(inputVector,outputVector,beta,dimension);
}

// for std vector in place
template <class ScalarType> void ComputeHouseholderVector(std::vector <ScalarType> &vector, double &beta)
{
    ComputeHouseholderVector(vector,vector,beta);
}

/**********************************/

template <class VectorType> double ComputeAngle(const VectorType &v1, const VectorType &v2)
{
    double normV1 = ComputeNorm(v1);
    double normV2 = ComputeNorm(v2);

    double t = ComputeScalarProduct(v1, v2) / (normV1 * normV2);

    if (t < -1.0)
        t = -1.0;
    if (t > 1.0)
        t = 1.0;

    return std::acos(t) * 180.0 / M_PI;
}

template <class VectorType> double ComputeOrientationAngle(const VectorType &v1, const VectorType &v2)
{
    double normV1 = ComputeNorm(v1);
    double normV2 = ComputeNorm(v2);

    double t = std::abs(ComputeScalarProduct(v1, v2)) / (normV1 * normV2);

    if (t > 1.0)
        t = 1.0;

    return std::acos(t) * 180.0 / M_PI;
}

template <class VectorType> void Revert(const VectorType &v, const unsigned int vSize, VectorType &resVec)
{
    for (unsigned int i = 0;i < vSize;++i)
        resVec[i] = -v[i];
}

template <class ScalarType> void Revert(const std::vector <ScalarType> &v, std::vector<ScalarType> &resVec)
{
    unsigned int length = v.size();
    Revert(v, length, resVec);
}

template <class ScalarType> void Revert(const itk::Point <ScalarType> &v, itk::Point <ScalarType> &resVec)
{
    unsigned int length = v.GetPointDimension();
    Revert(v, length, resVec);
}

template <class ScalarType, unsigned int NDimension>
void Revert(const vnl_vector_fixed <ScalarType, NDimension> &v, vnl_vector_fixed <ScalarType, NDimension> &resVec)
{
    Revert(v, NDimension, resVec);
}

template <class VectorType> void TransformCartesianToSphericalCoordinates(const VectorType &v, VectorType &resVec)
{
    double normV = ComputeNorm(v);
    Normalize(v, resVec);

    if (resVec[2] >= 1.0)
    {
        resVec[0] = 0;
        resVec[1] = 0;
    }
    else if (resVec[2] <= -1.0)
    {
        resVec[0] = M_PI;
        resVec[1] = 0;
    }
    else
    {
        double theta = std::acos(resVec[2]);
        double phi = std::atan2(resVec[1],resVec[0]);

        if (phi < 0.0)
            phi += (2.0 * M_PI);

        resVec[0] = theta;
        resVec[1] = phi;
    }

    resVec[2] = normV;
}

template <class VectorType> void TransformSphericalToCartesianCoordinates(const VectorType &v, VectorType &resVec)
{
    anima::TransformSphericalToCartesianCoordinates(v[0],v[1],v[2],resVec);
}

template <class VectorType> void TransformSphericalToCartesianCoordinates(double theta, double phi, double vectorNorm, VectorType &resVec)
{
    double tmpVal = std::cos(theta);
    if (std::abs(tmpVal) < 1e-4)
    {
        resVec[0] = std::cos(phi);
        resVec[1] = std::sin(phi);
        resVec[2] = 0;
    }
    else
    {
        resVec[0] = std::sin(theta) * std::cos(phi);
        resVec[1] = std::sin(theta) * std::sin(phi);
        resVec[2] = tmpVal;
    }

    for (unsigned int i = 0;i < 3;++i)
        resVec[i] *= vectorNorm;
}

template <class T1, class T2, class T3> void ProjectOnOrthogonalPlane(const std::vector <T1> &v, const std::vector <T2> &normalVec, std::vector <T3> &resVec)
{
    unsigned int length = v.size();

    if (normalVec.size() != length)
        throw itk::ExceptionObject(__FILE__, __LINE__,"Projection on orthogonal plane requires a normal vector of the same size as the input vector.",ITK_LOCATION);

    resVec.resize(length);

    double t = ComputeScalarProduct(normalVec, v);

    for (unsigned int i = 0;i < length;++i)
        resVec[i] = v[i] - t * normalVec[i];

    Normalize(resVec, resVec);
}

template <class T1, class T2, class T3, class T4> void RotateAroundAxis(const std::vector <T1> &v, const T2 &phi, const std::vector <T3> &normalVec, std::vector <T4> &resVec)
{
    unsigned int length = v.size();

    if (length != 3)
        throw itk::ExceptionObject(__FILE__, __LINE__,"Rotation around axis requires 3D input vector.",ITK_LOCATION);

    if (normalVec.size() != length)
        throw itk::ExceptionObject(__FILE__, __LINE__,"Rotation around axis requires a normal vector of the same size as the input vector.",ITK_LOCATION);

    resVec.resize(3);
    fill(resVec.begin(),resVec.end(),0);

    vnl_matrix <double> A(3,3);
    A.fill(0);
    A(0,1) = -normalVec[2];
    A(0,2) = normalVec[1];
    A(1,2) = -normalVec[0];
    A -= A.transpose();

    vnl_matrix <double> N(3,3);
    for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 3; ++j)
            N(i,j) = normalVec[i] * normalVec[j];

    vnl_matrix <double> R(3,3);
    R.fill(0.0);
    R.fill_diagonal(cos(phi));

    R += N * (1.0 - cos(phi)) + A * sin(phi);

    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0; j < 3; ++j)
            resVec[i] += R(i,j) * v[j];
}

template <class Vector3DType, class ScalarType> void RotateAroundAxis(const Vector3DType &v, const ScalarType &phi, const Vector3DType &normalVec, Vector3DType &resVec)
{
    vnl_matrix_fixed <ScalarType,3,3> A;
    A.fill(0.0);
    A(0,1) = -normalVec[2];
    A(0,2) = normalVec[1];
    A(1,2) = -normalVec[0];
    A -= A.transpose();

    vnl_matrix_fixed <ScalarType,3,3> N;
    for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 3; ++j)
            N(i,j) = normalVec[i] * normalVec[j];

    vnl_matrix_fixed <ScalarType,3,3> R;
    R.fill(0.0);
    R.fill_diagonal(std::cos(phi));

    R += N * (1.0 - std::cos(phi)) + A * std::sin(phi);

    for (unsigned int i = 0;i < 3;++i)
    {
        resVec[i] = 0;
        for (unsigned int j = 0;j < 3;++j)
            resVec[i] += R(i, j) * v[j];
    }
}

} // end of namespace anima
