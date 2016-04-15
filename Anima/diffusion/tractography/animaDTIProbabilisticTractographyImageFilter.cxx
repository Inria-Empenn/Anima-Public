#include "animaDTIProbabilisticTractographyImageFilter.h"
#include <cmath>

#include <animaVectorOperations.h>
#include <animaDistributionSampling.h>
#include <animaVMFDistribution.h>
#include <animaWatsonDistribution.h>
#include <animaLogarithmFunctions.h>

#include <itkSymmetricEigenAnalysis.h>

namespace anima
{

DTIProbabilisticTractographyImageFilter::DTIProbabilisticTractographyImageFilter()
: BaseProbabilisticTractographyImageFilter()
{
    m_ThresholdForProlateTensor = 0.25;
    m_OblateSigma = 0.1;
    
    m_FAThreshold = 0.2;
    
    this->SetModelDimension(6);
}

DTIProbabilisticTractographyImageFilter::~DTIProbabilisticTractographyImageFilter()
{
}

void DTIProbabilisticTractographyImageFilter::PrepareTractography()
{
    // Call base preparation
    BaseProbabilisticTractographyImageFilter::PrepareTractography();
    
    this->SetDesignMatrix();
    this->SetEstimationMatrix();
}

DTIProbabilisticTractographyImageFilter::Vector3DType
DTIProbabilisticTractographyImageFilter::ProposeNewDirection(Vector3DType &oldDirection, VectorType &modelValue,
                                                             Vector3DType &sampling_direction, double &log_prior,
                                                             double &log_proposal, std::mt19937 &random_generator,
                                                             unsigned int threadId)
{
    Vector3DType resVec(0.0);
    
    double concentrationParameter;
    bool is2d = this->GetInputImage(0)->GetLargestPossibleRegion().GetSize()[2] <= 1;
    double LC = this->GetLinearCoefficient(modelValue);
    
    if (LC > m_ThresholdForProlateTensor)
    {

        this->GetDTIPrincipalDirection(modelValue, sampling_direction, is2d);
                
        if (anima::ComputeScalarProduct(oldDirection, sampling_direction) < 0)
            sampling_direction *= -1;
        
        double FA = this->GetFractionalAnisotropy(modelValue);
                
        concentrationParameter = this->GetKappaFromFA(FA);
    }
    else
    {
        sampling_direction = oldDirection;
        concentrationParameter = this->GetKappaOfPriorDistribution();
    }
    
//    if (concentrationParameter > 700)
//        anima::SampleFromVMFDistributionNumericallyStable(concentrationParameter,sampling_direction,resVec,random_generator);
//    else
//        anima::SampleFromVMFDistribution(concentrationParameter,sampling_direction,resVec,random_generator);

    anima::SampleFromWatsonDistribution(concentrationParameter,sampling_direction,resVec,3,random_generator);
    
    if (is2d)
    {
        resVec[InputImageType::ImageDimension - 1] = 0;
        resVec.Normalize();
    }
    
    if (LC > m_ThresholdForProlateTensor)
    {
//        log_prior = anima::safe_log(anima::ComputeVMFPdf(resVec, oldDirection, this->GetKappaOfPriorDistribution()));
        log_prior = anima::safe_log(anima::EvaluateWatsonPDF(resVec, oldDirection, this->GetKappaOfPriorDistribution()));

//        log_proposal = anima::safe_log(anima::ComputeVMFPdf(resVec, sampling_direction, concentrationParameter));
        log_proposal = anima::safe_log(anima::EvaluateWatsonPDF(resVec, sampling_direction, concentrationParameter));
    }
    
    if (anima::ComputeScalarProduct(oldDirection, resVec) < 0)
        resVec *= -1;

    return resVec;
}

DTIProbabilisticTractographyImageFilter::Vector3DType
DTIProbabilisticTractographyImageFilter::
InitializeFirstIterationFromModel(Vector3DType &colinearDir, VectorType &modelValue, unsigned int threadId)
{    
    Vector3DType resVec, tmpVec;
    bool is2d = (this->GetInputImage(0)->GetLargestPossibleRegion().GetSize()[2] == 1);
	
    itk::SymmetricEigenAnalysis <Matrix3DType,Vector3DType,Matrix3DType> EigenAnalysis(InputImageType::ImageDimension);
    EigenAnalysis.SetOrderEigenValues(true);

    Matrix3DType dtiTensor, eigVecs;
    Vector3DType eigVals;

    unsigned int pos = 0;
    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j <= i;++j)
        {
            dtiTensor(i,j) = modelValue[pos];
            if (j != i)
                dtiTensor(j,i) = dtiTensor(i,j);
            ++pos;
        }

    EigenAnalysis.ComputeEigenValuesAndVectors(dtiTensor,eigVals,eigVecs);

    switch(this->GetInitialDirectionMode())
    {
        case Colinear:
        {
            double maxVal = 0;
            for (unsigned int i = 0;i < 3;++i)
            {
                for (unsigned int j = 0;j < 3;++j)
                    tmpVec[j] = eigVecs(i,j);

                double tmpVal = anima::ComputeScalarProduct(colinearDir,tmpVec);

                if (tmpVal < 0)
                {
                    tmpVec *= -1;
                    tmpVal *= -1;
                }

                if (maxVal < tmpVal)
                {
                    maxVal = tmpVal;
                    resVec = tmpVec;
                }
            }
            break;
        }

        case Weight:
        default:
        {
            for (unsigned int i = 0;i < 3;++i)
                resVec[i] = eigVecs(2,i);

            if (anima::ComputeScalarProduct(colinearDir, resVec) < 0)
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

bool DTIProbabilisticTractographyImageFilter::CheckModelProperties(double estimatedB0Value, double estimatedNoiseValue, VectorType &modelValue, unsigned int threadId)
{
    if (estimatedB0Value < 50.0)
        return false;
    
    bool isTensorNull = true;
    for (unsigned int j = 0;j < this->GetModelDimension();++j)
    {
        if (modelValue[j] != 0)
        {
            isTensorNull = false;
            break;
        }
    }
    
    if (isTensorNull)
        return false;

    double fractionalAnisotropy = this->GetFractionalAnisotropy(modelValue);
    if (fractionalAnisotropy < m_FAThreshold)
        return false;
    
    return true;
}

double DTIProbabilisticTractographyImageFilter::ComputeLogWeightUpdate(double b0Value, double noiseValue, Vector3DType &newDirection, Vector3DType &sampling_direction,
                                                                       VectorType &modelValue, VectorType &dwiValue,
                                                                       double &log_prior, double &log_proposal, unsigned int threadId)
{
    double resVal = 0;
    
    unsigned int numInputs = this->GetDiffusionGradients().size();
    bool is2d = this->GetInputImage(0)->GetLargestPossibleRegion().GetSize()[2] <= 1;
        
    // Computes prior, proposal and log-likelihood values
    double logLikelihood = 0;
    double LC = this->GetLinearCoefficient(modelValue);
    
    if (LC > m_ThresholdForProlateTensor)
    {
        Vector3DType dtiPrincipalDirection(0.0);
        
        this->GetDTIPrincipalDirection(modelValue, dtiPrincipalDirection, is2d);
        
        // Computes scalar values from DTI
        double perpLambda, meanLambda;
        this->GetEigenValueCombinations(modelValue,meanLambda,perpLambda);

        if (anima::ComputeScalarProduct(newDirection, dtiPrincipalDirection) < 0)
        	dtiPrincipalDirection *= -1;
        
        for (unsigned int i = 0;i < numInputs;++i)
        {
			double diffSignalValue = 0, meanSignalValue = 0;
            double var = noiseValue;

			if (var > 0)
			{
                double sclPrd = anima::ComputeScalarProduct(this->GetDiffusionGradient(i), dtiPrincipalDirection);
                double signalValue = b0Value * exp(-this->GetBValueItem(i) * (perpLambda + 3.0 * sclPrd * sclPrd * (meanLambda - perpLambda)));
                sclPrd = anima::ComputeScalarProduct(this->GetDiffusionGradient(i), newDirection);
                double rotatedSignalValue = b0Value * exp(-this->GetBValueItem(i) * (perpLambda + 3.0 * sclPrd * sclPrd * (meanLambda - perpLambda)));
	        
	            diffSignalValue = rotatedSignalValue - signalValue;
	            meanSignalValue = (rotatedSignalValue + signalValue) / 2.0;
	        
	            logLikelihood += (diffSignalValue * (dwiValue[i] - meanSignalValue) / var);
	    	}
        }
    }
    else
    {
        double localSigma = m_OblateSigma * m_OblateSigma;
        
        Vector3DType dtiMinorDirection(0.0);
        
        this->GetDTIMinorDirection(modelValue, dtiMinorDirection);
        
        double tmpVal = anima::ComputeScalarProduct(newDirection, dtiMinorDirection);
        if (tmpVal > 1.0)
            tmpVal = 1.0;
        if (tmpVal < -1.0)
            tmpVal = -1.0;
        
        tmpVal = std::acos(tmpVal) - M_PI / 2.0;
        
        logLikelihood -= tmpVal * tmpVal / (2.0 * localSigma);
    }
    
    resVal += log_prior - log_proposal;
    resVal += logLikelihood;
    
//    std::cout << "weight update: " << resVal << std::endl;
    
    return resVal;
}

void DTIProbabilisticTractographyImageFilter::SetDesignMatrix()
{
    unsigned int numInputs = this->GetDiffusionGradients().size();
    unsigned int numParameters = this->GetModelDimension() + 1;
    
    m_DesignMatrix.set_size(numInputs, numParameters);
    m_DesignMatrix.fill(0.0);
    m_DesignMatrix(0,0) = 1;
    
    for (unsigned int i = 0;i < numInputs;++i)
    {
        m_DesignMatrix(i,0) = 1;
        
        unsigned int pos = 1;
        for (unsigned int j = 0;j < 3;++j)
            for (unsigned int k = 0;k <= j;++k)
            {
                if (j != k)
                    m_DesignMatrix(i,pos) = - 2 * this->GetBValueItem(i) * this->GetDiffusionGradient(i)[j] * this->GetDiffusionGradient(i)[k];
                else
                    m_DesignMatrix(i,pos) = - this->GetBValueItem(i) * this->GetDiffusionGradient(i)[j] * this->GetDiffusionGradient(i)[j];
                
                ++pos;
            }
    }
}

void DTIProbabilisticTractographyImageFilter::SetEstimationMatrix()
{
    MatrixFreeType transposedDesignMatrix = m_DesignMatrix.transpose();
    MatrixFreeType squaredDesignMatrix = transposedDesignMatrix * m_DesignMatrix;
    
    m_EstimationMatrix = vnl_matrix_inverse<double>(squaredDesignMatrix).inverse() * transposedDesignMatrix;
}

void DTIProbabilisticTractographyImageFilter::SetKappaPolynomialCoefficients(std::vector <double> &coefs)
{
    m_KappaPolynomialCoefficients.resize(coefs.size());
    for (unsigned int i = 0;i < coefs.size();++i)
        m_KappaPolynomialCoefficients[i] = coefs[i];
}

double DTIProbabilisticTractographyImageFilter::GetKappaFromFA(double FA)
{
    double resVal = 0.0;
    for (unsigned int i = 0;i < m_KappaPolynomialCoefficients.size();++i)
        resVal += m_KappaPolynomialCoefficients[i] * std::pow(FA, (double)i);
    
    return 0.5 * resVal;
}

void DTIProbabilisticTractographyImageFilter::GetDTIPrincipalDirection(const VectorType &modelValue, Vector3DType &resVec, bool is2d)
{
    itk::SymmetricEigenAnalysis <Matrix3DType,Vector3DType,Matrix3DType> EigenAnalysis(InputImageType::ImageDimension);
    EigenAnalysis.SetOrderEigenValues(true);
    
    Matrix3DType dtiTensor, eigVecs;
    Vector3DType eigVals;
    
    unsigned int pos = 0;
    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j <= i;++j)
        {
            dtiTensor(i,j) = modelValue[pos];
            if (j != i)
                dtiTensor(j,i) = dtiTensor(i,j);
            ++pos;
        }
    
    EigenAnalysis.ComputeEigenValuesAndVectors(dtiTensor,eigVals,eigVecs);
    
    for (unsigned int i = 0;i < 3;++i)
    	resVec[i] = eigVecs(2,i);
    
    if (is2d)
    {
        resVec[2] = 0;
        resVec.Normalize();
    }
}
 
void DTIProbabilisticTractographyImageFilter::GetDTIMinorDirection(VectorType &modelValue, Vector3DType &resVec)
{
    itk::SymmetricEigenAnalysis <Matrix3DType,Vector3DType,Matrix3DType> EigenAnalysis(InputImageType::ImageDimension);
    EigenAnalysis.SetOrderEigenValues(true);
    
    Matrix3DType dtiTensor, eigVecs;
    Vector3DType eigVals;
    
    unsigned int pos = 0;
    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j <= i;++j)
        {
            dtiTensor(i,j) = modelValue[pos];
            if (j != i)
                dtiTensor(j,i) = dtiTensor(i,j);
            ++pos;
        }
    
    EigenAnalysis.ComputeEigenValuesAndVectors(dtiTensor,eigVals,eigVecs);
    
    for (unsigned int i = 0;i < 3;++i)
    	resVec[i] = eigVecs(0,i);
}

double DTIProbabilisticTractographyImageFilter::ComputeModelEstimation(DWIInterpolatorPointerVectorType &dwiInterpolators, ContinuousIndexType &index,
                                                                       VectorType &dwiValue, double &noiseValue, VectorType &modelValue)
{
	unsigned int numInputs = dwiInterpolators.size();
	
	dwiValue.SetSize(numInputs);
	
	for (unsigned int i = 0;i < numInputs;++i)
	    dwiValue[i] = dwiInterpolators[i]->EvaluateAtContinuousIndex(index);
    
    modelValue.SetSize(this->GetModelDimension());
    
    // Hard coded noise value for now
    noiseValue = 20;

    std::vector <double> logDWISignal(numInputs);
    
    // Search for the minimal b0 value
    double b0Value = 1.0e10;
    for (unsigned int i = 0;i < numInputs;++i)
    {
        if (this->GetBValueItem(i) == 0)
        {
            double tmpVal = dwiValue[i];
            
            if (tmpVal <= 0)
                tmpVal = 1.0e-4;
            
            logDWISignal[i] = anima::safe_log(tmpVal);
            
            if (tmpVal < b0Value)
                b0Value = tmpVal;
        }
    }
    
    // Compute log-signals
    for (unsigned int i = 0;i < numInputs;++i)
    {
        if (this->GetBValueItem(i) > 0)
        {
            double tmpVal = dwiValue[i];
            
            if (tmpVal > b0Value)
                tmpVal = b0Value - 1.0e-4;
            
            if (tmpVal <= 0)
                tmpVal = 1.0e-4;
            
            logDWISignal[i] = anima::safe_log(tmpVal);
        }
    }
    
    for (unsigned int i = 0;i < this->GetModelDimension();++i)
    {
        modelValue[i] = 0;
        for (unsigned int j = 0;j < numInputs;++j)
            modelValue[i] += m_EstimationMatrix(i+1,j) * logDWISignal[j];
    }
    
    double returnValue = 0;
    for (unsigned int j = 0;j < numInputs;++j)
        returnValue += m_EstimationMatrix(0,j) * logDWISignal[j];
    
    return exp(returnValue);
}

double DTIProbabilisticTractographyImageFilter::GetLinearCoefficient(VectorType &modelValue)
{
    itk::SymmetricEigenAnalysis <Matrix3DType,Vector3DType,Matrix3DType> EigenAnalysis(InputImageType::ImageDimension);
    EigenAnalysis.SetOrderEigenValues(true);
    
    Matrix3DType dtiTensor;
    Vector3DType eigVals;
    
    unsigned int pos = 0;
    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j <= i;++j)
        {
            dtiTensor(i,j) = modelValue[pos];
            if (j != i)
                dtiTensor(j,i) = dtiTensor(i,j);
            ++pos;
        }
    
    EigenAnalysis.ComputeEigenValues(dtiTensor,eigVals);
    
    double denom = 0;
    for (unsigned int i = 0;i < 3;++i)
        denom += eigVals[i] * eigVals[i];
    
    return (eigVals[2] - eigVals[1]) / sqrt(denom);
}

double DTIProbabilisticTractographyImageFilter::GetFractionalAnisotropy(VectorType &modelValue)
{
    itk::SymmetricEigenAnalysis <Matrix3DType,Vector3DType,Matrix3DType> EigenAnalysis(InputImageType::ImageDimension);
    EigenAnalysis.SetOrderEigenValues(true);
    
    Matrix3DType dtiTensor;
    Vector3DType eigVals;
    
    unsigned int pos = 0;
    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j <= i;++j)
        {
            dtiTensor(i,j) = modelValue[pos];
            if (j != i)
                dtiTensor(j,i) = dtiTensor(i,j);
            ++pos;
        }
    
    EigenAnalysis.ComputeEigenValues(dtiTensor,eigVals);
    
    double meanLambda = 0;
    for (unsigned int i = 0;i < 3;++i)
        meanLambda += eigVals[i];

    meanLambda /= 3.0;

    double num = 0;
    double denom = 0;
    for (unsigned int i = 0;i < 3;++i)
    {
        num += (eigVals[i] - meanLambda) * (eigVals[i] - meanLambda);
        denom += eigVals[i] * eigVals[i];
    }
    
    // FA
    return sqrt(3.0 * num / (2.0 * denom));
}

void DTIProbabilisticTractographyImageFilter::GetEigenValueCombinations(VectorType &modelValue, double &meanLambda, double &perpLambda)
{
    itk::SymmetricEigenAnalysis <Matrix3DType,Vector3DType,Matrix3DType> EigenAnalysis(InputImageType::ImageDimension);
    EigenAnalysis.SetOrderEigenValues(true);
    
    Matrix3DType dtiTensor;
    Vector3DType eigVals;
    
    unsigned int pos = 0;
    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j <= i;++j)
        {
            dtiTensor(i,j) = modelValue[pos];
            if (j != i)
                dtiTensor(j,i) = dtiTensor(i,j);
            ++pos;
        }
    
    EigenAnalysis.ComputeEigenValues(dtiTensor,eigVals);
    
    meanLambda = 0;
    perpLambda = 0;
    for (unsigned int i = 0;i < 3;++i)
    {
        meanLambda += eigVals[i];
        if (i == 1)
            perpLambda = meanLambda / 2.0;
    }
    
    meanLambda /= 3.0;
}

} // end of namespace anima
