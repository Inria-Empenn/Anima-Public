#include "animaODFProbabilisticTractographyImageFilter.h"
#include <cmath>
#include <random>

#include <animaODFMaximaCostFunction.h>
#include <animaNewuoaOptimizer.h>

#include <animaVectorOperations.h>
#include <animaMatrixOperations.h>
#include <animaDistributionSampling.h>
#include <animaVMFDistribution.h>
#include <animaWatsonDistribution.h>

namespace anima
{

ODFProbabilisticTractographyImageFilter::ODFProbabilisticTractographyImageFilter()
: BaseProbabilisticTractographyImageFilter()
{
    m_ODFSHOrder = 4;
    m_GFAThreshold = 0.1;
    m_MinimalDiffusionProbability = 0;

    // Magic value for regularization from Descoteaux et al. article
    m_Lambda = 0.006;

    // Delta from Aganj et al article to better handle noise
    m_DeltaAganjRegularization = 0.001;

    m_CurvatureScale = 6.0;

    m_ODFSHBasis = NULL;

    this->SetModelDimension(15);
}

ODFProbabilisticTractographyImageFilter::~ODFProbabilisticTractographyImageFilter()
{
    if (m_ODFSHBasis)
        delete m_ODFSHBasis;
}

void ODFProbabilisticTractographyImageFilter::SetODFSHOrder(unsigned int num)
{
    m_ODFSHOrder = num;
    this->SetModelDimension((m_ODFSHOrder + 1)*(m_ODFSHOrder + 2)/2);
}

void ODFProbabilisticTractographyImageFilter::PrepareTractography()
{
    // Call base preparation
    BaseProbabilisticTractographyImageFilter::PrepareTractography();

    // Initialize estimation matrices for Aganj et al based estimation
    if (m_ODFSHBasis)
        delete m_ODFSHBasis;

    m_ODFSHBasis = new anima::ODFSphericalHarmonicBasis(m_ODFSHOrder);

    unsigned int posValue = 0;
    unsigned int numGrads = this->GetDiffusionGradients().size();

    vnl_matrix <double> BMatrix(numGrads,this->GetModelDimension());
    Vector3DType sphDiffGradient;

    unsigned int pos = 0;
    for (unsigned int i = 0;i < numGrads;++i)
    {
        if (anima::ComputeNorm(this->GetDiffusionGradient(i)) == 0)
            continue;

        posValue = 0;
        anima::TransformCartesianToSphericalCoordinates(this->GetDiffusionGradient(i),sphDiffGradient);
        for (int k = 0;k <= (int)m_ODFSHOrder;k += 2)
            for (int m = -k;m <= k;++m)
            {
                BMatrix(pos,posValue) = m_ODFSHBasis->getNthSHValueAtPosition(k,m,sphDiffGradient[0],sphDiffGradient[1]);
                ++posValue;
            }

        ++pos;
    }

    std::vector <double> LVector(this->GetModelDimension(),0);
    std::vector <double> PVector(this->GetModelDimension(),0);

    posValue = 0;
    for (unsigned int k = 0;k <= m_ODFSHOrder;k += 2)
    {
        double ljVal = k*k*(k + 1)*(k + 1);
        double pjValNum = 1, pjValDenom = 1, pjVal;

        for (unsigned int l = 2;l <= k; l += 2)
        {
            pjValNum *= (l-1);
            pjValDenom *= l;
        }

        pjVal = pjValNum/pjValDenom;
        if (k/2 % 2 != 0)
            pjVal *= -1;

        for (int m = -k;m <= (int)k;++m)
        {
            LVector[posValue] = ljVal;
            PVector[posValue] = k*(k+1.0)*pjVal/(-8.0*M_PI);

            ++posValue;
        }
    }

    vnl_matrix <double> tmpMat = BMatrix.transpose() * BMatrix;
    for (unsigned int i = 0;i < this->GetModelDimension();++i)
        tmpMat(i,i) += m_Lambda*LVector[i];

    vnl_matrix_inverse <double> tmpInv(tmpMat);

    m_SignalCoefsMatrix = tmpInv.inverse() * BMatrix.transpose();
    unsigned int realNumGrads = pos;
    m_TMatrix.set_size(this->GetModelDimension(),realNumGrads);

    for (unsigned int i = 0;i < this->GetModelDimension();++i)
        for (unsigned int j = 0;j < realNumGrads;++j)
            m_TMatrix(i,j) = m_SignalCoefsMatrix(i,j) * PVector[i];
}

ODFProbabilisticTractographyImageFilter::Vector3DType
ODFProbabilisticTractographyImageFilter::ProposeNewDirection(Vector3DType &oldDirection, VectorType &modelValue,
                                                             Vector3DType &sampling_direction, double &log_prior,
                                                             double &log_proposal, std::mt19937 &random_generator,
                                                             unsigned int threadId)
{
    Vector3DType resVec(0.0);
    bool is2d = (this->GetInputImage(0)->GetLargestPossibleRegion().GetSize()[2] == 1);

    DirectionVectorType maximaODF;
    unsigned int numDirs = this->FindODFMaxima(modelValue,maximaODF,m_MinimalDiffusionProbability,is2d);
    ListType mixtureWeights(numDirs,0);
    ListType kappaValues(numDirs,0);

    double chosenKappa = 0;

    Vector3DType sphDirection;

    double sumWeights = 0;

    for (unsigned int i = 0;i < numDirs;++i)
    {
        if (anima::ComputeScalarProduct(oldDirection, maximaODF[i]) < 0)
            maximaODF[i] *= -1;

        anima::TransformCartesianToSphericalCoordinates(maximaODF[i],sphDirection);
        mixtureWeights[i] = m_ODFSHBasis->getValueAtPosition(modelValue,sphDirection[0],sphDirection[1]);

        // 0.5 is for Watson kappa
        kappaValues[i] = 0.5 * m_CurvatureScale * m_ODFSHBasis->getCurvatureAtPosition(modelValue,sphDirection[0],sphDirection[1]);

        if ((std::isnan(kappaValues[i]))||(kappaValues[i] <= 0)||(kappaValues[i] >= 1000))
            mixtureWeights[i] = 0;

        sumWeights += mixtureWeights[i];
    }

    if (sumWeights == 0)
    {
        numDirs = 0;
        sampling_direction = oldDirection;
        chosenKappa = this->GetKappaOfPriorDistribution();
    }
    else
    {
        for (unsigned int i = 0;i < numDirs;++i)
            mixtureWeights[i] /= sumWeights;

        std::discrete_distribution<> dist(mixtureWeights.begin(),mixtureWeights.end());
        unsigned int chosenDirection = dist(random_generator);

        sampling_direction = maximaODF[chosenDirection];
        chosenKappa = kappaValues[chosenDirection];
    }

//    if (chosenKappa > 700)
//        anima::SampleFromVMFDistributionNumericallyStable(chosenKappa,sampling_direction,resVec,random_generator);
//    else
//        anima::SampleFromVMFDistribution(chosenKappa,sampling_direction,resVec,random_generator);

    anima::SampleFromWatsonDistribution(chosenKappa,sampling_direction,resVec,3,random_generator);

    if (is2d)
    {
        resVec[InputImageType::ImageDimension - 1] = 0;
        resVec.Normalize();
    }

    if (numDirs > 0)
    {
//        log_prior = anima::safe_log( anima::ComputeVMFPdf(resVec, oldDirection, this->GetKappaOfPriorDistribution()));
        log_prior = anima::safe_log( anima::EvaluateWatsonPDF(resVec, oldDirection, this->GetKappaOfPriorDistribution()));

        log_proposal = 0;
        for (unsigned int i = 0;i < numDirs;++i)
        {
//            log_proposal += mixtureWeights[i] * anima::ComputeVMFPdf(resVec, maximaODF[i], kappaValues[i]);
            log_proposal += mixtureWeights[i] * anima::EvaluateWatsonPDF(resVec, maximaODF[i], kappaValues[i]);
        }

        log_proposal = anima::safe_log(log_proposal);
    }

    if (anima::ComputeScalarProduct(oldDirection, resVec) < 0)
        resVec *= -1;

    return resVec;
}

ODFProbabilisticTractographyImageFilter::Vector3DType ODFProbabilisticTractographyImageFilter::InitializeFirstIterationFromModel(Vector3DType &colinearDir, VectorType &modelValue,
                                                                                                                                 unsigned int threadId)
{
    Vector3DType resVec(0.0), tmpVec;
    bool is2d = (this->GetInputImage(0)->GetLargestPossibleRegion().GetSize()[2] == 1);

    DirectionVectorType maximaODF;
    unsigned int numDirs = this->FindODFMaxima(modelValue,maximaODF,m_MinimalDiffusionProbability,is2d);
    if (numDirs == 0)
        return colinearDir;

    switch (this->GetInitialDirectionMode())
    {
        case Colinear:
        {
            double maxVal = 0;
            for (unsigned int i = 0;i < numDirs;++i)
            {
                for (unsigned int j = 0;j < InputImageType::ImageDimension;++j)
                    tmpVec[j] = maximaODF[i][j];

                double tmpVal = anima::ComputeScalarProduct(colinearDir, tmpVec);

                if (tmpVal < 0)
                {
                    tmpVec *= -1;
                    tmpVal *= -1;
                }

                if (tmpVal > maxVal)
                {
                    resVec = tmpVec;
                    maxVal = tmpVal;
                }
            }

            break;
        }

        case Weight:
        default:
        {
            resVec = maximaODF[0];
            if (anima::ComputeScalarProduct(colinearDir,resVec) < 0)
                resVec *= -1;
            break;
        }
    }

    if (is2d)
    {
        resVec[2] = 0;
        resVec.Normalize();
    }

    return resVec;
}

bool ODFProbabilisticTractographyImageFilter::CheckModelProperties(double estimatedB0Value, double estimatedNoiseValue, VectorType &modelValue, unsigned int threadId)
{
    if (estimatedB0Value < 50.0)
        return false;

    bool isModelNull = true;
    for (unsigned int j = 0;j < this->GetModelDimension();++j)
    {
        if (modelValue[j] != 0)
        {
            isModelNull = false;
            break;
        }
    }

    if (isModelNull)
        return false;

    double fractionalAnisotropy = this->GetGeneralizedFractionalAnisotropy(modelValue);
    if (fractionalAnisotropy < m_GFAThreshold)
        return false;

    return true;
}

double ODFProbabilisticTractographyImageFilter::ComputeLogWeightUpdate(double b0Value, double noiseValue, Vector3DType &newDirection, Vector3DType &sampling_direction,
                                                                       VectorType &modelValue, VectorType &dwiValue,
                                                                       double &log_prior, double &log_proposal, unsigned int threadId)
{
    Matrix3DType rotationMatrix = anima::GetRotationMatrixFromVectors(sampling_direction,newDirection);
    rotationMatrix = rotationMatrix.GetTranspose();

    unsigned int numGrads = this->GetDiffusionGradients().size();
    Vector3DType rotatedGradient(0.0);
    Vector3DType sphRotatedGradient, tmpGradDir(0.0);
    std::vector <double> rescaledModel(this->GetModelDimension(),0);

    // Compute scale model to get signal simulation
    unsigned int realNumGrads = m_TMatrix.columns();
    std::vector <double> tmpData(realNumGrads,0);

    unsigned int pos = 0;
    for (unsigned int i = 0;i < numGrads;++i)
    {
        if (anima::ComputeNorm(this->GetDiffusionGradient(i)) == 0)
            continue;

        double e = dwiValue[i] / b0Value;

        if (e < 0)
            tmpData[pos] = m_DeltaAganjRegularization / 2.0;
        else if (e < m_DeltaAganjRegularization)
            tmpData[pos] = m_DeltaAganjRegularization / 2.0 + e * e / (2.0 * m_DeltaAganjRegularization);
        else if (e < 1.0 - m_DeltaAganjRegularization)
            tmpData[pos] = e;
        else if (e < 1)
            tmpData[pos] = 1.0 - m_DeltaAganjRegularization / 2.0 - (1.0 - e) * (1.0 - e) / (2.0 * m_DeltaAganjRegularization);
        else
            tmpData[pos] = 1.0 - m_DeltaAganjRegularization / 2.0;

        tmpData[pos] = std::log(-std::log(tmpData[pos]));
        ++pos;
    }

    for (unsigned int i = 0;i < this->GetModelDimension();++i)
    {
        for (unsigned int j = 0;j < realNumGrads;++j)
            rescaledModel[i] += m_SignalCoefsMatrix(i,j)*tmpData[j];
    }

    double logLikelihood = 0;

    for (unsigned int i = 0;i < numGrads;++i)
    {
        if (anima::ComputeNorm(this->GetDiffusionGradient(i)) == 0)
            continue;

        for (unsigned int j = 0;j < InputImageType::ImageDimension;++j)
        {
            rotatedGradient[j] = 0;
            for (unsigned int k = 0;k < InputImageType::ImageDimension;++k)
                rotatedGradient[j] += rotationMatrix(j,k) * this->GetDiffusionGradient(i)[k];
        }

        anima::TransformCartesianToSphericalCoordinates(rotatedGradient,sphRotatedGradient);

        double signalValue = m_ODFSHBasis->getValueAtPosition(rescaledModel,sphRotatedGradient[0],sphRotatedGradient[1]);

        // Aganj et al formulation
        signalValue = b0Value * exp(-exp(signalValue));
        //std::cout << signalValue << " " << dwiValue[i] << std::endl;

        double tmpVal = - 0.5 * log(noiseValue * 2 * M_PI);
        double residual = dwiValue[i] - signalValue;
        tmpVal -= residual * residual / (2.0 * noiseValue);

        logLikelihood += tmpVal;
    }

    double resVal = logLikelihood / numGrads;
    resVal += log_prior - log_proposal;

    return resVal;
}

//! Returns ODF maxima ordered by probability of diffusion, only those with probability superior to minVal are given
unsigned int ODFProbabilisticTractographyImageFilter::FindODFMaxima(const VectorType &modelValue, DirectionVectorType &maxima, double minVal, bool is2d)
{
    ListType modelValueList(modelValue.GetSize());
    for (unsigned int i = 0;i < modelValue.GetSize();++i)
        modelValueList[i] = modelValue[i];

    std::vector < std::vector <double> > initDirs(15);

    // Find the max of the ODF v, set max to the value...
    struct XYZ array[] = {
        {-1,0,0},
        {0,-1,0},
        {0,0,-1},
        {0.1789,0.1113,-0.9776},
        {0.0635,-0.3767,-0.9242},
        {-0.7108,-0.0516,-0.7015},
        {-0.6191,0.4385,-0.6515},
        {-0.2424,-0.7843,-0.571},
        {0.2589,0.618,-0.7423},
        {0.8169,-0.1697,-0.5513},
        {0.8438,-0.5261,-0.106},
        {0.2626,-0.9548,-0.1389},
        {-1e-04,-0.9689,0.2476},
        {-0.7453,-0.6663,0.0242},
        {-0.9726,-0.2317,0.0209}
    };

    for (unsigned int i = 0;i < initDirs.size();++i)
    {
        initDirs[i].resize(3);
        initDirs[i][0] = array[i].x;
        initDirs[i][1] = array[i].y;
        initDirs[i][2] = array[i].z;
    }

    typedef anima::ODFMaximaCostFunction CostFunctionType;
    typedef anima::NewuoaOptimizer OptimizerType;

    OptimizerType::Pointer opt = OptimizerType::New();
    typedef std::map < double, Vector3DType > MapType;
    MapType dmap;

    std::vector <double> angle;
    OptimizerType::ParametersType tmpValue(2);
    Vector3DType tmpValueVector(1.0);
    Vector3DType cartesianVector(1.0);

    for (unsigned int i = 0; i < initDirs.size();++i)
    {
        anima::TransformCartesianToSphericalCoordinates(initDirs[i],angle);
        tmpValue[0] = angle[0];
        tmpValue[1] = angle[1];

        CostFunctionType::Pointer cost = CostFunctionType::New();
        cost->SetODFSHOrder(m_ODFSHOrder);
        cost->SetBasisParameters(modelValueList);

        opt->SetCostFunction(cost);

        opt->SetRhoBegin(M_PI/80);
        opt->SetRhoEnd(M_PI/200);
        opt->SetNumberSamplingPoints(cost->GetNumberOfParameters() + 2);
        opt->SetMaximumIteration(100);
        opt->SetMaximize(true);

        opt->SetInitialPosition(tmpValue);
        opt->StartOptimization();

        tmpValue = opt->GetCurrentPosition();
        tmpValueVector[0] = tmpValue[0];
        tmpValueVector[1] = tmpValue[1];

        // Some check needed here to see if we really found a maximum
        anima::TransformSphericalToCartesianCoordinates(tmpValueVector,cartesianVector);
        dmap[opt->GetValue()] = cartesianVector;
    }

    // Find true maximas
    std::vector <bool> usefulMaxima(15,true);
    unsigned int pos = 0;
    for (MapType::reverse_iterator it = dmap.rbegin();it != dmap.rend();++it)
    {
        unsigned int posIn = pos+1;
        MapType::reverse_iterator itin = it;
        ++itin;
        for (;itin != dmap.rend();++itin)
        {
            if (anima::ComputeOrientationAngle((*it).second,(*itin).second) < 15)
                usefulMaxima[posIn] = false;

            ++posIn;
        }

        ++pos;
    }

    pos = 0;
    maxima.clear();
    for (MapType::reverse_iterator it = dmap.rbegin();it != dmap.rend();++it)
    {
        if ((usefulMaxima[pos])&&((*it).first > minVal))
            maxima.push_back((*it).second);

        ++pos;
    }

    if (is2d)
    {
        std::vector <bool> outOfPlaneDirs(maxima.size(),false);
        for (unsigned int i = 0;i < maxima.size();++i)
        {
            maxima[i][2] = 0;

            double norm = 0;
            for (unsigned int j = 0;j < InputImageType::ImageDimension - 1;++j)
                norm += maxima[i][j] * maxima[i][j];
            norm = sqrt(norm);

            outOfPlaneDirs[i] = (fabs(norm) < 0.5);

            if (!outOfPlaneDirs[i])
            {
                for (unsigned int j = 0;j < InputImageType::ImageDimension - 1;++j)
                    maxima[i][j] /= norm;
            }
        }

        DirectionVectorType outMaxima;

        for (unsigned int i = 0;i < maxima.size();++i)
        {
            if (!outOfPlaneDirs[i])
                outMaxima.push_back(maxima[i]);
        }

        maxima = outMaxima;
    }

    return maxima.size();
}

double ODFProbabilisticTractographyImageFilter::ComputeModelEstimation(DWIInterpolatorPointerVectorType &dwiInterpolators, ContinuousIndexType &index,
                                                                       VectorType &dwiValue, double &noiseValue, VectorType &modelValue)
{
    unsigned int numInputs = dwiInterpolators.size();
    dwiValue.SetSize(numInputs);

    for (unsigned int i = 0;i < numInputs;++i)
        dwiValue[i] = dwiInterpolators[i]->EvaluateAtContinuousIndex(index);

    double b0Value = dwiValue[0];

    modelValue.SetSize(this->GetModelDimension());
    modelValue.Fill(0.0);

    // Hard coded noise value for now
    noiseValue = 20;
    if (b0Value <= 0)
        return 0;

    unsigned int realNumGrads = m_TMatrix.columns();
    std::vector <double> tmpData(realNumGrads,0);

    unsigned int pos = 0;
    for (unsigned int i = 0;i < numInputs;++i)
    {
        if (anima::ComputeNorm(this->GetDiffusionGradient(i)) == 0)
            continue;

        double e = dwiValue[i] / b0Value;

        if (e < 0)
            tmpData[pos] = m_DeltaAganjRegularization / 2.0;
        else if (e < m_DeltaAganjRegularization)
            tmpData[pos] = m_DeltaAganjRegularization / 2.0 + e * e / (2.0 * m_DeltaAganjRegularization);
        else if (e < 1.0 - m_DeltaAganjRegularization)
            tmpData[pos] = e;
        else if (e < 1)
            tmpData[pos] = 1.0 - m_DeltaAganjRegularization / 2.0 - (1.0 - e) * (1.0 - e) / (2.0 * m_DeltaAganjRegularization);
        else
            tmpData[pos] = 1.0 - m_DeltaAganjRegularization / 2.0;

        tmpData[pos] = std::log(-std::log(tmpData[pos]));
        ++pos;
    }

    modelValue[0] = 1.0 / (2.0 * sqrt(M_PI));

    for (unsigned int i = 1;i < this->GetModelDimension();++i)
    {
        for (unsigned int j = 0;j < realNumGrads;++j)
            modelValue[i] += m_TMatrix(i,j)*tmpData[j];
    }

    return b0Value;
}

double ODFProbabilisticTractographyImageFilter::GetGeneralizedFractionalAnisotropy(VectorType &modelValue)
{
    double sumSquares = 0;
    for (unsigned int i = 0;i < this->GetModelDimension();++i)
        sumSquares += modelValue[i]*modelValue[i];

    return sqrt(1.0 - modelValue[0]*modelValue[0]/sumSquares);
}

} // end of namespace anima
