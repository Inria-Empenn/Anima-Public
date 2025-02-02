#pragma once
#include <animaQuaternionTools.h>

namespace anima
{

template <class TInput, class TOutput> vnl_matrix <TOutput> computeRotationFromQuaternion(vnl_vector <TInput> eigenVector)
{
    double normVector = 0;

    for (unsigned int i = 0;i < eigenVector.size();++i)
        normVector += eigenVector[i]*eigenVector[i];

    vnl_matrix <TOutput> resVal(3,3,0);

    resVal(0,0) = (1.0/normVector) * (eigenVector[0]*eigenVector[0] + eigenVector[1]*eigenVector[1] - eigenVector[2]*eigenVector[2] - eigenVector[3]*eigenVector[3]);
    resVal(0,1) = (2.0/normVector) * (eigenVector[1]*eigenVector[2] - eigenVector[0]*eigenVector[3]);
    resVal(0,2) = (2.0/normVector) * (eigenVector[1]*eigenVector[3] + eigenVector[0]*eigenVector[2]);

    resVal(1,0) = (2.0/normVector) * (eigenVector[1]*eigenVector[2] + eigenVector[0]*eigenVector[3]);
    resVal(1,1) = (1.0/normVector) * (eigenVector[0]*eigenVector[0] - eigenVector[1]*eigenVector[1] + eigenVector[2]*eigenVector[2] - eigenVector[3]*eigenVector[3]);
    resVal(1,2) = (2.0/normVector) * (eigenVector[2]*eigenVector[3] - eigenVector[0]*eigenVector[1]);

    resVal(2,0) = (2.0/normVector) * (eigenVector[1]*eigenVector[3] - eigenVector[0]*eigenVector[2]);
    resVal(2,1) = (2.0/normVector) * (eigenVector[2]*eigenVector[3] + eigenVector[0]*eigenVector[1]);
    resVal(2,2) = (1.0/normVector) * (eigenVector[0]*eigenVector[0] - eigenVector[1]*eigenVector[1] - eigenVector[2]*eigenVector[2] + eigenVector[3]*eigenVector[3]);

    return resVal;
}

template <class TInput, class TOutput, unsigned int NDimensions>
void pairingToQuaternion(const vnl_vector_fixed <TInput,NDimensions> &inputPoint, const vnl_vector_fixed <TInput,NDimensions> &inputTransformedPoint,
                         vnl_matrix <TOutput> &outputMatrix)
{
    return pairingToQuaternion(inputPoint,inputTransformedPoint,outputMatrix,NDimensions);
}

template <class TInput, class TOutput, unsigned int NDimensions>
void pairingToQuaternion(const itk::Vector <TInput,NDimensions> &inputPoint, const itk::Vector <TInput,NDimensions> &inputTransformedPoint,
                         vnl_matrix <TOutput> &outputMatrix)
{
    return pairingToQuaternion(inputPoint,inputTransformedPoint,outputMatrix,NDimensions);
}

template <class PointType, class TOutput>
void pairingToQuaternion(const PointType &inputPoint, const PointType &inputTransformedPoint, vnl_matrix <TOutput> &outputMatrix, unsigned int ndim)
{
    outputMatrix.set_size(ndim+1,ndim+1);
    outputMatrix.fill(0);

    for (unsigned int i = 0;i < ndim;++i)
    {
        switch (i)
        {
            case 0:
                outputMatrix(0,1) = inputPoint[i] - inputTransformedPoint[i];
                outputMatrix(2,3) = - (inputPoint[i] + inputTransformedPoint[i]);
                outputMatrix(1,0) = - outputMatrix(0,1);
                outputMatrix(3,2) = - outputMatrix(2,3);
                break;

            case 1:
                outputMatrix(0,2) = inputPoint[i] - inputTransformedPoint[i];
                outputMatrix(1,3) = inputPoint[i] + inputTransformedPoint[i];
                outputMatrix(2,0) = - outputMatrix(0,2);
                outputMatrix(3,1) = - outputMatrix(1,3);
                break;

            case 2:
                outputMatrix(0,3) = inputPoint[i] - inputTransformedPoint[i];
                outputMatrix(1,2) = - (inputPoint[i] + inputTransformedPoint[i]);
                outputMatrix(2,1) = - outputMatrix(1,2);
                outputMatrix(3,0) = - outputMatrix(0,3);
                break;

            default:
                break;
        }
    }
}

template <class TInput, class TOutput, unsigned int NDimensions>
void pairingToQuaternionScalingsDerivative(const vnl_vector_fixed <TInput, NDimensions> &inputPoint, const vnl_vector_fixed <TInput, NDimensions> &inputTransformedPoint,
                                           vnl_matrix <TOutput> &outputMatrix, const int &dimScal)
{
    return pairingToQuaternionScalingsDerivative(inputPoint, inputTransformedPoint, outputMatrix, NDimensions, dimScal);
}

template <class TInput, class TOutput, unsigned int NDimensions>
void pairingToQuaternionScalingsDerivative(const itk::Vector <TInput, NDimensions> &inputPoint, const itk::Vector <TInput, NDimensions> &inputTransformedPoint,
                                           vnl_matrix <TOutput> &outputMatrix, const int &dimScal)
{
    return pairingToQuaternionScalingsDerivative(inputPoint, inputTransformedPoint, outputMatrix, NDimensions, dimScal);
}

template <class PointType, class TOutput>
void pairingToQuaternionScalingsDerivative(const PointType &inputPoint, const PointType &inputTransformedPoint, vnl_matrix <TOutput> &outputMatrix, unsigned int ndim, const int &dimScal)
{
    outputMatrix.set_size(ndim + 1, ndim + 1);
    outputMatrix.fill(0);

    switch (dimScal)
    {
        case 0:
            outputMatrix(0, 2) = -inputTransformedPoint[2];
            outputMatrix(0, 3) = inputTransformedPoint[1];
            outputMatrix(1, 2) = inputTransformedPoint[1];
            outputMatrix(1, 3) = inputTransformedPoint[2];
            outputMatrix += outputMatrix.transpose();
            outputMatrix(0, 0) = inputTransformedPoint[0];
            outputMatrix(1, 1) = inputTransformedPoint[0];
            outputMatrix(2, 2) = -inputTransformedPoint[0];
            outputMatrix(3, 3) = -inputTransformedPoint[0];
            outputMatrix *= inputPoint[0];
            break;

        case 1:
            outputMatrix(0, 1) = inputTransformedPoint[2];
            outputMatrix(0, 3) = -inputTransformedPoint[0];
            outputMatrix(1, 2) = inputTransformedPoint[0];
            outputMatrix(2, 3) = inputTransformedPoint[2];
            outputMatrix += outputMatrix.transpose();
            outputMatrix(0, 0) = inputTransformedPoint[1];
            outputMatrix(1, 1) = -inputTransformedPoint[1];
            outputMatrix(2, 2) = inputTransformedPoint[1];
            outputMatrix(3, 3) = -inputTransformedPoint[1];
            outputMatrix *= inputPoint[1];
            break;

        case 2:
            outputMatrix(0, 1) = -inputTransformedPoint[1];
            outputMatrix(0, 2) = inputTransformedPoint[0];
            outputMatrix(1, 3) = inputTransformedPoint[0];
            outputMatrix(2, 3) = inputTransformedPoint[1];
            outputMatrix += outputMatrix.transpose();
            outputMatrix(0, 0) = inputTransformedPoint[2];
            outputMatrix(1, 1) = -inputTransformedPoint[2];
            outputMatrix(2, 2) = -inputTransformedPoint[2];
            outputMatrix(3, 3) = inputTransformedPoint[2];
            outputMatrix *= inputPoint[2];
            break;

        default:
            break;
    }
    
}

} // end of namespace anima
