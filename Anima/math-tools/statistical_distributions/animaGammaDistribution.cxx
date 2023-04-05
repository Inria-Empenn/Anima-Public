#include "animaGammaDistribution.h"
#include <cmath>
#include <limits>
#include <itkMacro.h>

namespace anima
{

double GammaDistribution::GetDensity(const double &x)
{
    return std::exp(this->GetLogDensity(x));
}

double GammaDistribution::GetLogDensity(const double &x)
{
    double shapeParameter = this->GetShapeParameter();
    double scaleParameter = this->GetScaleParameter();
    double resValue = (shapeParameter - 1.0) * std::log(x);
    resValue -= x / scaleParameter;
    resValue -= std::lgamma(shapeParameter);
    resValue -= shapeParameter * std::log(scaleParameter);
    return resValue;
}

void GammaDistribution::Fit(const VectorType &sample, const std::string &method)
{
    unsigned int dimValue = sample.size();
    double shapeParameter, scaleParameter;

    if (method == "mle")
    {
        // Code for estimating theta and kappa via MLE
        double meanValue = 0.0;
        double meanLogValue = 0.0;
        for (unsigned int i = 0;i < dimValue;++i)
        {
            double tmpValue = sample[i];
            meanValue += tmpValue;
            meanLogValue += std::log(tmpValue);
        }
        meanValue /= (double)dimValue;
        meanLogValue /= (double)dimValue;
        double logMeanValue = std::log(meanValue);

        double sValue = logMeanValue - meanLogValue;
        shapeParameter = (3.0 - sValue + std::sqrt((sValue - 3.0) * (sValue - 3.0) + 24.0 * sValue)) / (12.0 * sValue);
        scaleParameter = 0.0;
        if (shapeParameter > std::numeric_limits<double>::epsilon())
            scaleParameter = meanValue / shapeParameter;
    }
    else if (method == "biased-closed-form" || method == "unbiased-closed-form")
    {
        // Code for estimating theta and kappa via CF biased and unbiased
        double sumValue = 0.0;
        double sumLogValue = 0.0;
        double sumLogXValue = 0.0;
        for (unsigned int i = 0;i < dimValue;++i)
        {
            double tmpValue = sample[i];
            sumValue += tmpValue;
            sumLogValue += std::log(tmpValue);
            sumLogXValue += (std::log(tmpValue) * tmpValue);
        }

        double denValue = (double)dimValue * sumLogXValue - sumLogValue * sumValue;
        shapeParameter = (double)dimValue * sumValue;
        scaleParameter = (double)dimValue * sumLogXValue - sumLogValue * sumValue;
        if (denValue > std::numeric_limits<double>::epsilon())
            shapeParameter /= denValue;
            scaleParameter /= ((double)dimValue * (double)dimValue);
        
        if (method == "unbiased-closed-form")
        {
            shapeParameter -= (3.0 * shapeParameter - 2.0 / 3.0 * shapeParameter / (1.0 + shapeParameter) - 4.0 / 5.0 * shapeParameter / ((1.0 + shapeParameter) * (1.0 + shapeParameter))) / (double)dimValue;
            scaleParameter *= (double)dimValue / ((double)dimValue - 1.0);
        }
    }
    else
        throw itk::ExceptionObject(__FILE__, __LINE__, "Unsupported estimation method for the parameters of the Gamma distribution.", ITK_LOCATION);

    this->SetShapeParameter(shapeParameter);
    this->SetScaleParameter(scaleParameter);
}

void GammaDistribution::Random(VectorType &sample, GeneratorType &generator)
{
    DistributionType distributionValue(this->GetShapeParameter(), this->GetScaleParameter());
    unsigned int nSamples = sample.size();
    for (unsigned int i = 0;i < nSamples;++i)
        sample[i] = distributionValue(generator);
}
    
} // end of namespace anima
