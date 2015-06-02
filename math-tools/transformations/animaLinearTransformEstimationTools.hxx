#pragma once
#include "animaLinearTransformEstimationTools.h"

#include <animaMatrixLogExp.h>

#include <vnl/vnl_inverse.h>
#include <itkVector.h>

#include <itkMatrixOffsetTransformBase.h>
#include <itkTranslationTransform.h>
#include <itkContinuousIndex.h>

namespace anima
{
template <class TInput, class TScalarType, unsigned int NDimensions>
void computeTranslationLSWFromTranslations(std::vector < itk::Point<TInput,NDimensions> > &inputOrigins,
                                           std::vector < itk::Point<TInput,NDimensions> > &inputTransformed,
                                           std::vector <TInput> &weights,
                                           typename itk::AffineTransform<TScalarType,NDimensions>::Pointer &resultTransform)
{
    unsigned int nbPts = inputOrigins.size();
    itk::Point <TInput,NDimensions> barX, barY;
    // See Pennec PhD
    vnl_matrix <TInput> AMatrix(4,4,0);
    vnl_matrix <TInput> tmpMatrix(4,4);

    for (unsigned int j = 0;j < NDimensions;++j)
    {
        barX[j] = 0;
        barY[j] = 0;
    }

    double sumWeights = 0;
    for (unsigned int i = 0;i < nbPts;++i)
    {
        for (unsigned int j = 0;j < NDimensions;++j)
        {
            barX[j] += weights[i] * inputOrigins[i][j];
            barY[j] += weights[i] * inputTransformed[i][j];
        }

        sumWeights += weights[i];
    }

    for (unsigned int j = 0;j < NDimensions;++j)
    {
        barX[j] /= sumWeights;
        barY[j] /= sumWeights;
    }

    for (unsigned int i = 0;i < nbPts;++i)
    {
        anima::pairingToQuaternion(inputOrigins[i] - barX,inputTransformed[i] - barY,tmpMatrix);
        AMatrix += weights[i]*tmpMatrix.transpose()*tmpMatrix;
    }

    itk::SymmetricEigenAnalysis < vnl_matrix <TInput>, vnl_diag_matrix<TInput>, vnl_matrix <TInput> > eigenSystem(4);
    vnl_matrix <double> eVec(4,4);
    vnl_diag_matrix <double> eVals(4);

    eigenSystem.SetOrderEigenValues(true);
    eigenSystem.ComputeEigenValuesAndVectors(AMatrix, eVals, eVec);

    vnl_matrix <TScalarType> rotationMatrix = anima::computeRotationFromQuaternion<TInput,TScalarType>(eVec.get_row(0));

    itk::Vector <TScalarType,NDimensions> translationPart;
    for (unsigned int i = 0;i < NDimensions;++i)
        translationPart[i] = barY[i] - barX[i];

    resultTransform = itk::AffineTransform<TScalarType,NDimensions>::New();

    resultTransform->SetIdentity();
    resultTransform->SetOffset(translationPart);
}

template <class TInput, class TScalarType, unsigned int NDimensions>
void computeRigidLSWFromTranslations(std::vector < itk::Point<TInput,NDimensions> > &inputOrigins,
                                     std::vector < itk::Point<TInput,NDimensions> > &inputTransformed,
                                     std::vector <TInput> &weights,
                                     typename itk::AffineTransform<TScalarType,NDimensions>::Pointer &resultTransform)
{
    unsigned int nbPts = inputOrigins.size();
    itk::Point <TInput,NDimensions> barX, barY;
    // See Pennec PhD
    vnl_matrix <TInput> AMatrix(4,4,0);
    vnl_matrix <TInput> tmpMatrix(4,4);

    for (unsigned int j = 0;j < NDimensions;++j)
    {
        barX[j] = 0;
        barY[j] = 0;
    }

    double sumWeights = 0;
    for (unsigned int i = 0;i < nbPts;++i)
    {
        for (unsigned int j = 0;j < NDimensions;++j)
        {
            barX[j] += weights[i] * inputOrigins[i][j];
            barY[j] += weights[i] * inputTransformed[i][j];
        }

        sumWeights += weights[i];
    }

    for (unsigned int j = 0;j < NDimensions;++j)
    {
        barX[j] /= sumWeights;
        barY[j] /= sumWeights;
    }

    for (unsigned int i = 0;i < nbPts;++i)
    {
        anima::pairingToQuaternion(inputOrigins[i] - barX,inputTransformed[i] - barY,tmpMatrix);
        AMatrix += weights[i]*tmpMatrix.transpose()*tmpMatrix;
    }

    itk::SymmetricEigenAnalysis < vnl_matrix <TInput>, vnl_diag_matrix<TInput>, vnl_matrix <TInput> > eigenSystem(4);
    vnl_matrix <double> eVec(4,4);
    vnl_diag_matrix <double> eVals(4);

    eigenSystem.SetOrderEigenValues(true);
    eigenSystem.ComputeEigenValuesAndVectors(AMatrix, eVals, eVec);

    vnl_matrix <TScalarType> rotationMatrix = anima::computeRotationFromQuaternion<TInput,TScalarType>(eVec.get_row(0));

    itk::Vector <TScalarType,NDimensions> translationPart;
    for (unsigned int i = 0;i < NDimensions;++i)
    {
        translationPart[i] = barY[i];
        for (unsigned int j = 0;j < NDimensions;++j)
            translationPart[i] -= rotationMatrix(i,j)*barX[j];
    }

    resultTransform = itk::AffineTransform<TScalarType,NDimensions>::New();

    resultTransform->SetMatrix(rotationMatrix);
    resultTransform->SetOffset(translationPart);
}

template <class TInput, class TScalarType, unsigned int NDimensions>
void computeLogEuclideanAverage(std::vector < vnl_matrix <TInput> > &inputTransforms, std::vector <TInput> &weights,
                                typename itk::AffineTransform<TScalarType,NDimensions>::Pointer &resultTransform)
{
    unsigned int nbPts = inputTransforms.size();

    vnl_matrix <TInput> resultMatrix(NDimensions+1,NDimensions+1,0);
    vnl_matrix <TInput> tmpMatrix;
    double sumWeights = 0;
    for (unsigned int i = 0;i < nbPts;++i)
    {
        sumWeights += weights[i];
        resultMatrix += weights[i] * inputTransforms[i];
    }

    resultMatrix /= sumWeights;

    resultMatrix = anima::GetExponential(resultMatrix);

    resultTransform = itk::AffineTransform<TScalarType,NDimensions>::New();

    vnl_matrix <TScalarType> affinePart(NDimensions,NDimensions);
    itk::Vector <TScalarType,NDimensions> translationPart;

    for (unsigned int i = 0;i < NDimensions;++i)
    {
        translationPart[i] = resultMatrix(i,NDimensions);
        for (unsigned int j = 0;j < NDimensions;++j)
            affinePart(i,j) = resultMatrix(i,j);
    }

    resultTransform->SetMatrix(affinePart);
    resultTransform->SetOffset(translationPart);
}

template <class TInput, class TScalarType, unsigned int NDimensions>
void computeAffineLSWFromTranslations(std::vector < itk::Point<TInput,NDimensions> > &inputOrigins,
                                      std::vector < itk::Point<TInput,NDimensions> > &inputTransformed,
                                      std::vector <TInput> &weights,
                                      typename itk::AffineTransform<TScalarType,NDimensions>::Pointer &resultTransform)
{
    unsigned int nbPts = inputOrigins.size();
    itk::Point <TInput,NDimensions> barX, barY;
    // See Pennec PhD
    vnl_matrix <TInput> SigmaXXMatrix(NDimensions,NDimensions,0);
    vnl_matrix <TInput> SigmaYXMatrix(NDimensions,NDimensions,0);

    for (unsigned int j = 0;j < NDimensions;++j)
    {
        barX[j] = 0;
        barY[j] = 0;
    }

    double sumWeights = 0;
    for (unsigned int i = 0;i < nbPts;++i)
    {
        for (unsigned int j = 0;j < NDimensions;++j)
        {
            barX[j] += weights[i] * inputOrigins[i][j];
            barY[j] += weights[i] * inputTransformed[i][j];
        }

        sumWeights += weights[i];
    }

    for (unsigned int j = 0;j < NDimensions;++j)
    {
        barX[j] /= sumWeights;
        barY[j] /= sumWeights;
    }

    for (unsigned int i = 0;i < nbPts;++i)
    {
        for (unsigned int j = 0;j < NDimensions;++j)
            for (unsigned int k = 0;k < NDimensions;++k)
            {
                SigmaXXMatrix(j,k) += weights[i]*(inputOrigins[i][j] - barX[j])*(inputOrigins[i][k] - barX[k]);
                SigmaYXMatrix(j,k) += weights[i]*(inputTransformed[i][j] - barY[j])*(inputOrigins[i][k] - barX[k]);
            }
    }

    vnl_matrix <TInput> affineMatrix = SigmaYXMatrix * vnl_inverse (SigmaXXMatrix);

    vnl_matrix <TScalarType> outMatrix(NDimensions,NDimensions);
    itk::Vector <TScalarType,NDimensions> translationPart;
    for (unsigned int i = 0;i < NDimensions;++i)
    {
        translationPart[i] = barY[i];
        for (unsigned int j = 0;j < NDimensions;++j)
        {
            outMatrix(i,j) = affineMatrix(i,j);
            translationPart[i] -= affineMatrix(i,j)*barX[j];
        }
    }

    resultTransform = itk::AffineTransform<TScalarType,NDimensions>::New();

    resultTransform->SetMatrix(outMatrix);
    resultTransform->SetOffset(translationPart);
}

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

template <class TInput, class TOutput, unsigned int NDimensions> void pairingToQuaternion(itk::Vector<TInput,NDimensions> inputPoint,
                                                                                          itk::Vector<TInput,NDimensions> inputTransformedPoint,
                                                                                          vnl_matrix <TOutput> &outputMatrix)
{
    outputMatrix.set_size(4,4);
    outputMatrix.fill(0);

    for (unsigned int i = 0;i < NDimensions;++i)
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

} // end of namespace anima
