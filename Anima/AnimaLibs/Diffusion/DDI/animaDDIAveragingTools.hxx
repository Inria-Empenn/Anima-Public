#pragma once
#include "animaDDIAveragingTools.h"

#include <animaMatrixLogExp.h>
#include <animaMatrixOperations.h>
#include <animaBaseTensorTools.h>

#include <itkSymmetricEigenAnalysis.h>

namespace anima
{

template <class ScalarType> void ComputeVoxelWeights(const typename itk::VectorImage<ScalarType, 3>::IndexType &inputIndex,
                                                     std::vector<typename itk::VectorImage<ScalarType, 3>::IndexType> &variableIndex,
                                                     std::vector<ScalarType> &wPosition, const typename itk::Image<ScalarType, 3>::SizeType &bound,
                                                     const int step, unsigned int nbVoxels)
{
    int boundx = bound[0];
    int boundy = bound[1];
    int boundz = bound[2];

    int xIndex = inputIndex[0];
    int yIndex = inputIndex[1];
    int zIndex = inputIndex[2];

    int modx = xIndex % step;
    int mody = yIndex % step;
    int modz = zIndex % step;

    double scale = 1/double(step);

    typename itk::Image<ScalarType, 3>::IndexType tempIndex;

    variableIndex.clear();
    wPosition.clear();

    //Enumerate all the cases
    if (xIndex < boundx - step + 1 && yIndex < boundy - step + 1 && zIndex < boundz - step + 1)
    {
        // Voxel on the grid
        if (modx == 0 && mody == 0 && modz == 0)
            nbVoxels = 0;

        // Two voxels
        else if (modx != 0 && mody == 0 && modz == 0)
        {
            nbVoxels = 2;
            tempIndex[1] = yIndex;
            tempIndex[2] = zIndex;

            for (int varx = - modx; varx < step; varx += step)
            {
                tempIndex[0] = xIndex + varx;
                variableIndex.push_back(tempIndex);
                wPosition.push_back(1 - scale*std::abs(varx));
            }
        }
        else if (modx == 0 && mody != 0 && modz == 0)
        {
            nbVoxels = 2;
            tempIndex[0] = xIndex;
            tempIndex[2] = zIndex;
            for (int vary = - mody; vary < step; vary += step)
            {
                tempIndex[1] = yIndex + vary;
                variableIndex.push_back(tempIndex);
                wPosition.push_back(1 - scale*std::abs(vary));
            }
        }
        else if (modx == 0 && mody == 0 && modz != 0)
        {
            nbVoxels = 2;
            tempIndex[0] = xIndex;
            tempIndex[1] = yIndex;
            for (int varz = - modz; varz < step; varz += step)
            {
                tempIndex[2] = zIndex + varz;
                variableIndex.push_back(tempIndex);
                wPosition.push_back(1 - scale*std::abs(varz));
            }
        }

        // Four voxels
        else if (modx != 0 && mody != 0 && modz == 0)
        {
            nbVoxels = 4;
            tempIndex[2] = zIndex;
            for (int varx = -modx; varx < step; varx += step)
            {
                for (int vary = -mody; vary < step; vary += step)
                {
                    tempIndex[0] = xIndex + varx;
                    tempIndex[1] = yIndex + vary;
                    variableIndex.push_back(tempIndex);
                    wPosition.push_back((1 - scale*std::abs(varx))*(1 - scale*std::abs(vary)));
                }
            }
        }
        else if (modx != 0 && mody == 0 && modz != 0)
        {
            nbVoxels = 4;
            tempIndex[1] = yIndex;
            for (int varx = -modx; varx < step; varx += step)
            {
                for (int varz = -modz; varz < step; varz += step)
                {
                    tempIndex[0] = xIndex + varx;
                    tempIndex[2] = zIndex + varz;
                    variableIndex.push_back(tempIndex);
                    wPosition.push_back((1 - scale*std::abs(varx))*(1 - scale*std::abs(varz)));
                }
            }
        }
        else if (modx == 0 && mody != 0 && modz != 0)
        {
            nbVoxels = 4;
            tempIndex[0] = xIndex;
            for (int vary = -mody; vary < step; vary += step)
            {
                for (int varz = -modz; varz < step; varz += step)
                {
                    tempIndex[1] = yIndex + vary;
                    tempIndex[2] = zIndex + varz;
                    variableIndex.push_back(tempIndex);
                    wPosition.push_back((1 - scale*std::abs(vary))*(1 - scale*std::abs(varz)));
                }
            }
        }

        //Eight voxels
        else if (modx != 0 && mody != 0 && modz != 0)
        {
            nbVoxels = 8;
            for (int varx = -modx; varx < step; varx += step)
            {
                for (int vary = -mody; vary < step; vary += step)
                {
                    for (int varz = -modz; varz < step; varz += step)
                    {
                        tempIndex[0] = xIndex + modx;
                        tempIndex[1] = yIndex + mody;
                        tempIndex[2] = zIndex + modz;
                        variableIndex.push_back(tempIndex);
                        wPosition.push_back((1 - scale*std::abs(varx))*(1 - scale*std::abs(vary))*(1 - scale*std::abs(varz)));
                    }
                }
            }
        }
    }
    else
        nbVoxels = 0;

}

template<class ScalarType> void RemoveNullFascicles(std::vector<ScalarType>& v, std::vector<ScalarType>& d,
                                                    std::vector<ScalarType>& kappa,
                                                    std::vector<vnl_vector<ScalarType> >& directions,
                                                    std::vector<ScalarType>& w,
                                                    std::vector<ScalarType>& wFree,
                                                    const bool optionsBoundV)
{
    int nbFascicles = v.size();
    int nbFreeWater = wFree.size();

    double sumTotal = 0;

    for (int i = 0;i < nbFascicles;++i)
    {
        //If the weight is null remove the fascicule
        if (w[i] == 0)
        {
            v.erase(v.begin() + i);
            d.erase(d.begin() + i);
            kappa.erase(kappa.begin() + i);
            directions.erase(directions.begin() + i);
            w.erase(w.begin() + i);

            nbFascicles--;
            i--;
        }
        else
            sumTotal += w[i];
    }

    // Add free water area to the sumTotal
    for (int i = 0;i < nbFreeWater;++i)
        sumTotal += wFree[i];

    //Correct weights of fascicles
    for (int i = 0;i < nbFascicles;++i)
    {
        w[i] /= sumTotal;
        // Set boundaries of v between 0.01 and 0.99
        if (optionsBoundV)
        {
            if (v[i] == 1)
                v[i] = 0.99;
            if (v[i] == 0)
                v[i] = 0.01;
        }
    }

    //Correct weights of free water too.
    for (int i = 0; i < nbFreeWater; ++i)
        wFree[i] /= sumTotal;
}

template<class ScalarType> void CreateCovarianceMatrixFromDDIParameter(const std::vector<ScalarType>& v, const std::vector<ScalarType>& d,
                                                                       const std::vector<ScalarType>& kappa, const std::vector<vnl_vector<ScalarType> >& mu,
                                                                       std::vector<vnl_matrix<ScalarType> >& Sigma)

{
    int nbFascicles = v.size();

    vnl_matrix<ScalarType> U(3,3);
    vnl_matrix<ScalarType> Id(3,3);
    Id.set_identity();
    vnl_matrix<ScalarType> muMatrix(3,1);
    for (int i = 0; i < nbFascicles; ++i)
    {
        // Copy vector mu into a matrix (3,1) to do the product
        for (int j = 0; j < 3; j ++)
            muMatrix(j,0) = mu[i](j);

        U = muMatrix*muMatrix.transpose();
        Sigma[i] = (1 - v[i])*(d[i]/(1 + kappa[i])) * (Id + kappa[i]*U);

    }
}

template<class ScalarType> void DDIAveraging(std::vector<ScalarType>& v, std::vector<ScalarType>& d,
                                             std::vector<ScalarType>& kappa,
                                             std::vector<vnl_vector<ScalarType> >& directions,
                                             std::vector<ScalarType>& w,
                                             const int method, double &averageNu, double &averageDiffusivity,
                                             double &averageKappa, vnl_vector <ScalarType> &averageDirection,
                                             double &averageWeight)
{
    // Initialize some standard value
    ScalarType step = 1;
    ScalarType epsilon = 0.0001;

    // First of all test if there is some null input value and remove it.
    int nbFascicles = v.size();

    averageNu = 0;
    averageDiffusivity = 0;
    averageKappa = 0;
    averageDirection.fill(0);
    averageWeight = 0;

    if (nbFascicles == 0)
        return;

    if (nbFascicles == 1)
    {
        averageNu = v[0];
        averageDiffusivity = d[0];
        averageKappa = kappa[0];
        averageDirection = directions[0];
        averageWeight = w[0];

        return;
    }

    // Start by putting all directions into the same half sphere from direction cosine matrix
    vnl_matrix <double> dirsMatrix(3,3,0);
    vnl_matrix <double> tmpMat(3,3);
    double totalWeights = 0;

    anima::LogEuclideanTensorCalculator <double>::Pointer leCalculator = anima::LogEuclideanTensorCalculator <double>::New();

    for (unsigned int i = 0;i < nbFascicles;++i)
    {
        tmpMat.set_identity();
        tmpMat *= epsilon;

        for (unsigned int j = 0;j < 3;++j)
            for (unsigned int k = j;k < 3;++k)
            {
                double tmpVal = directions[i][j] * directions[i][k];
                tmpMat(j,k) += tmpVal;
                if (j != k)
                    tmpMat(k,j) += tmpVal;
            }

        leCalculator->GetTensorLogarithm(tmpMat,tmpMat);
        dirsMatrix += w[i] * tmpMat;
        totalWeights += w[i];
    }

    dirsMatrix /= totalWeights;
    leCalculator->GetTensorExponential(dirsMatrix,dirsMatrix);

    vnl_matrix <double> eigVecs(3,3);
    vnl_diag_matrix <double> eigVals(3);

    itk::SymmetricEigenAnalysis < vnl_matrix <double>, vnl_diag_matrix <double> > eigSystem(3);
    eigSystem.SetOrderEigenValues(true);
    eigSystem.ComputeEigenValuesAndVectors(dirsMatrix,eigVals,eigVecs);

    //Put all directions in the same half sphere
    for (unsigned int i = 0; i < nbFascicles; ++i)
    {
        double scalarProduct = directions[i][2];

        if (scalarProduct < 0)
            directions[i] *= -1;

        directions[i] /= directions[i].two_norm();
    }

    // Then create the covariance matrix from parameters
    std::vector<vnl_matrix<ScalarType> > Sigma (nbFascicles);
    CreateCovarianceMatrixFromDDIParameter(v, d, kappa, directions, Sigma);

    // Compute
    std::vector <ScalarType> radius;
    for (int i = 0; i< nbFascicles; i++)
        radius.push_back(std::sqrt(v[i]*d[i]));

    // Then normalize the sum of weight to one
    for (unsigned int i = 0;i < nbFascicles;++i)
        w[i] /= totalWeights;

    // Estimate directly radius
    ScalarType mRadius2 = 0;
    for (unsigned int i = 0;i < nbFascicles;++i)
        mRadius2 += w[i]*radius[i]*radius[i];

    //Estimate directly Sigma with the log-euclidean Frechet mean
    vnl_matrix <double> lSigma (3,3,0), mSigma(3,3);

    for (unsigned int i = 0;i < nbFascicles;++i)
    {
        leCalculator->GetTensorLogarithm(Sigma[i],tmpMat);
        lSigma += w[i]*tmpMat;
    }

    leCalculator->GetTensorExponential(lSigma,mSigma);

    //Get eigenvalue and eigen vector of mSigma
    vnl_matrix <double> sigmaEigVecs(3,3);
    vnl_diag_matrix <double> sigmaEigVals(3);
    eigSystem.ComputeEigenValuesAndVectors(mSigma,sigmaEigVals,sigmaEigVecs);

    ScalarType L1 = sigmaEigVals(2,2);
    ScalarType L2 = sigmaEigVals(1,1);
    ScalarType L3 = sigmaEigVals(0,0);
    ScalarType Lperp = std::sqrt(L2*L3);
    vnl_vector<ScalarType> muSigma = sigmaEigVecs.get_row(2);

    vnl_vector<ScalarType> mU(3,0);
    ScalarType mV1 = 0, mV2 = 0, mV = 0, mD = 0, mKappa = 0;

    //Estimate mu, kappa and nu following different methods
    switch (method)
    {
        //Simplest averaging (Classic)
        case 0:
        {
            std::vector<vnl_vector<ScalarType> > originalAngles(nbFascicles);
            vnl_vector<ScalarType> averageAngle(3, 0);
            averageAngle[0] = 1;
            for (int i = 0; i < nbFascicles; i++)
            {
                anima::TransformCartesianToSphericalCoordinates(directions[i], originalAngles[i]);
                averageAngle[1] += w[i]*originalAngles[i][1];
                averageAngle[2] += w[i]*originalAngles[i][2];

                mV += w[i]*v[i];
                mD += w[i]*d[i];
                mKappa += w[i]*kappa[i];
            }
            //Come back in cartesian coordinate
            anima::TransformSphericalToCartesianCoordinates(averageAngle, mU);
            break;
        }

        //Tensor averaging
        case 1:
        {
            mU = eigVecs.get_row(2);
            for (int i = 0; i < nbFascicles; i++)
            {
                mV += w[i]*v[i];
                mD += w[i]*d[i];
                mKappa += w[i]*kappa[i];
            }
            break;
        }

        //log VMF
        case 2:
        {
            double maxWeight = 0;
            double maxIndex = 0;
            for (unsigned int i = 0;i < nbFascicles;++i)
            {
                if (w[i] > maxWeight)
                {
                    maxWeight = w[i];
                    maxIndex = i;
                }
            }
            mKappa = kappa[maxIndex];
            ScalarType ConditionStop = 0.001;
            ScalarType lKappa = ConditionStop + 1;
            while (std::abs(lKappa) > ConditionStop )
            {
                lKappa = 0;
                //Average in the tangent space
                for (unsigned int i = 0; i < nbFascicles; ++i)
                    lKappa = lKappa + (step/nbFascicles)*w[i]*log(kappa[i] / mKappa);

                mKappa = mKappa*exp(lKappa);
            }
            mU = eigVecs.get_row(2);
            mV1 = mRadius2/(L1 + mRadius2);
            mV2 = mRadius2/(L2 * (mKappa + 1) + mRadius2);
            mV = (mV1 + mV2)/2;
            mD = mRadius2/mV;

            break;
        }

            //covariance Analytic
        case 3:
        {
            mU = muSigma;
            mV = mRadius2/(L1+mRadius2);
            mKappa = (L1 - Lperp)/Lperp;
            mD = mRadius2/mV;

            break;
        }
    }

    //Record results in the DDI format
    for (unsigned int i = 0;i < 3;++i)
        averageDirection[i] = mU[i];
    averageKappa = mKappa;
    averageDiffusivity = mD;
    averageNu = mV;
    averageWeight = totalWeights;
}

template <class ScalarType>
void
FreeWaterDDIAveraging(std::vector < std::vector <ScalarType> > &wFree,
                      std::vector < std::vector <ScalarType> > &dFree,
                      std::vector <ScalarType> &averageDiffusivities,
                      std::vector <ScalarType> &averageFreeWeights)

{
    unsigned int numIsoCompartments = wFree.size();
    unsigned int nbFascicles = 0;

    if (numIsoCompartments > 0)
        nbFascicles = wFree[0].size();

    averageDiffusivities.resize(numIsoCompartments);
    std::fill(averageDiffusivities.begin(),averageDiffusivities.end(),0.0);
    averageFreeWeights.resize(numIsoCompartments);
    std::fill(averageFreeWeights.begin(),averageFreeWeights.end(),0.0);

    for (unsigned int j = 0;j < numIsoCompartments;++j)
    {
        for (unsigned int i = 0; i < nbFascicles; ++i)
        {
            if ((wFree[j][i] != 0)&&(dFree[j][i] != 0))
            {
                averageFreeWeights[j] += wFree[j][i];
                averageDiffusivities[j] += wFree[j][i] * std::log(dFree[j][i]);
            }
        }
    }

    for (unsigned int j = 0;j < numIsoCompartments;++j)
    {
        if (averageFreeWeights[j] != 0)
            averageDiffusivities[j] = std::exp(averageDiffusivities[j] / averageFreeWeights[j]);
    }
}

template< class ScalarType>
void ComputeDistanceMatrixBetweenFascicles(std::vector<ScalarType>& v, std::vector<ScalarType>& d,
                                           std::vector<ScalarType>& kappa,
                                           std::vector<vnl_vector<ScalarType> >& directions,
                                           int method, vnl_matrix<ScalarType>& distanceMatrix)

{
    int kappaW = 1;
    int radiusW = 1;
    double epsilon = 0.0001;
    int nbFascicles = v.size();

    distanceMatrix.set_size(nbFascicles,nbFascicles);
    distanceMatrix.fill(0.0);

    std::vector <ScalarType> workV(2);
    std::vector <ScalarType> workD(2);
    std::vector <ScalarType> workKappa(2);
    std::vector <vnl_matrix<ScalarType> > workU(2);
    std::vector <vnl_vector<ScalarType> > workUVec(2);

    for (int i = 0; i < 2; ++i)
    {
        for (int k=0; k < 3; ++k)
        {
            workUVec[i].set_size(3);
            workU[i].set_size(3,1);
        }
    }
    std::vector < vnl_matrix<ScalarType> > workSigma(2);

    vnl_matrix <ScalarType> kappaMatrix(nbFascicles, nbFascicles, 0);
    vnl_matrix <ScalarType> radiusMatrix(nbFascicles, nbFascicles, 0);
    vnl_matrix <ScalarType> uMatrix(nbFascicles, nbFascicles, 0);

    double uSum = 0, kappaSum = 0, rayonSum = 0;
    distanceMatrix.fill(0);

    anima::LogEuclideanTensorCalculator <double>::Pointer leCalculator = anima::LogEuclideanTensorCalculator <double>::New();

    for (int i = 0; i < nbFascicles; ++i)
    {
        for (int j = i + 1; j < nbFascicles; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                workU[0][k][0] = directions[i][k];
                workU[1][k][0] = directions[j][k];
                workUVec[0][k] = directions[i][k];
                workUVec[1][k] = directions[j][k];
            }

            workKappa[0] = kappa[i];
            workKappa[1] = kappa[j];

            workD[0] = d[i];
            workD[1] = d[j];

            workV[0] = v[i];
            workV[1] = v[j];

            // Classic
            if (method == 0)
            {
                std::vector<vnl_vector<double> > angles(2);
                TransformCartesianToSphericalCoordinates(workUVec[0], angles[0]);
                TransformCartesianToSphericalCoordinates(workUVec[1], angles[1]);

                uMatrix[i][j] = std::sqrt(std::pow(angles[0][1] - angles[1][1], 2) + std::pow(angles[0][2] - angles[1][2], 2));
                radiusMatrix[i][j] = std::abs(workV[0]*workD[0] - workV[1]*workD[1]);
                kappaMatrix[i][j] = std::abs(workKappa[0] - workKappa[1]);

                uSum += uMatrix[i][j];
                kappaSum += kappaMatrix[i][j];
                rayonSum += radiusMatrix[i][j];

            }
            else if ((method == 1) || (method == 2)) // Tensor and log VMF
            {
                vnl_matrix<ScalarType> epsMatrix(3,3);
                epsMatrix.set_identity();
                epsMatrix *= epsilon;
                vnl_matrix<ScalarType> u0Tensor = epsMatrix + workU[0]*workU[0].transpose();
                vnl_matrix<ScalarType> u1Tensor = epsMatrix + workU[1]*workU[1].transpose();

                vnl_matrix<ScalarType> logu0(3,3);
                vnl_matrix<ScalarType> logu1(3,3);
                vnl_matrix<ScalarType> Diff(3,3,0);

                leCalculator->GetTensorLogarithm(u0Tensor,logu0);
                leCalculator->GetTensorLogarithm(u1Tensor,logu1);
                Diff = logu0 - logu1;
                uMatrix[i][j] = Diff.frobenius_norm();
                radiusMatrix[i][j] = std::abs(workV[0]*workD[0] - workV[1]*workD[1]);
                if (method == 2)
                    kappaMatrix[i][j] = std::abs(std::log(workKappa[0]) - std::log(workKappa[1]));
                else
                    kappaMatrix[i][j] = std::abs(workKappa[0] - workKappa[1]);

                uSum += uMatrix[i][j];
                kappaSum += kappaMatrix[i][j];
                rayonSum += radiusMatrix[i][j];
            }
            else if (method == 3) // covariance Analytic
            {
                CreateCovarianceMatrixFromDDIParameter(workV, workD, workKappa, workUVec, workSigma);
                vnl_matrix<ScalarType> logSigma0(3,3,0);
                vnl_matrix<ScalarType> logSigma1(3,3,0);
                vnl_matrix<ScalarType> Diff(3,3,0);

                leCalculator->GetTensorLogarithm(workSigma[0],logSigma0);
                leCalculator->GetTensorLogarithm(workSigma[1],logSigma1);
                Diff = logSigma0 - logSigma1;

                double distVal = Diff.frobenius_norm();
                distanceMatrix(i,j) = distVal * distVal;
            }
            else
                std::cerr << "Method needs to be included between 0 and 3" << std::endl;
        }
    }

    // Compute distance matrix for all cases exept analytic
    if (method != 3)
    {
        for (int i = 0; i < nbFascicles; ++i)
        {
            for (int j = i + 1; j < nbFascicles; ++j)
            {
                if (uSum != 0)
                    uMatrix[i][j] = uMatrix[i][j]/uSum;
                else
                    uMatrix[i][j] = 0;

                if (kappaSum != 0)
                    kappaMatrix[i][j] = kappaW*kappaMatrix[i][j]/kappaSum;
                else
                    kappaMatrix[i][j] = 0;

                if (rayonSum != 0)
                    radiusMatrix[i][j] = radiusW*radiusMatrix[i][j]/rayonSum;
                else
                    radiusMatrix[i][j] = 0;

                double distVal = uMatrix[i][j] + kappaMatrix[i][j] + radiusMatrix[i][j];
                distanceMatrix(i,j) = distVal * distVal;
            }
        }
    }

    for (int i = 0; i < nbFascicles; ++i)
    {
        for (int j = i + 1; j < nbFascicles; ++j)
            distanceMatrix(j,i) = distanceMatrix(i,j);
    }
}

} // end of namespace anima
