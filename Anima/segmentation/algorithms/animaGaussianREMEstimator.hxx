#include "animaGaussianREMEstimator.h"

namespace anima
{

template <typename TInputImage, typename TMaskImage>
bool GaussianREMEstimator<TInputImage,TMaskImage>::concentration()
{
    std::vector<vnl_matrix<GaussianREMEstimator::NumericType> > inverseCovariance;
    inverseCovariance.reserve(this->m_GaussianModel.size());
    double *determinantCovariance = new double[this->m_GaussianModel.size()];

    if(!this->m_APosterioriProbability.empty())
        this->m_APosterioriProbability.clear();
    if(!this->m_JointHistogram.empty())
        this->m_JointHistogram.clear();

    GaussianFunctionType::CovarianceMatrixType *covari = new GaussianFunctionType::CovarianceMatrixType[this->m_GaussianModel.size()];
    for(unsigned int i = 0; i < this->m_GaussianModel.size(); i++)
    {
        covari[i] = (this->m_GaussianModel[i])->GetCovariance();
    }

    //1. We calculate covariance inverse and determinant
    for(unsigned int i = 0; i < this->m_GaussianModel.size(); i++)
    {
        determinantCovariance[i] = vnl_determinant(covari[i].GetVnlMatrix());
        if(std::fabs(determinantCovariance[i]) < 1e-12)
        {
            return false;
        }
        inverseCovariance.push_back(covari[i].GetInverse());
    }

    Histogram::iterator histoIt;
    GaussianFunctionType::CovarianceMatrixType intensities(this->m_OriginalJointHistogram.begin()->first.size(),1);
    std::vector<double> probas(this->m_GaussianModel.size());

    ResidualMap residualMap;
    double numberOfPixels=0;

    GaussianFunctionType::CovarianceMatrixType *xi = new GaussianFunctionType::CovarianceMatrixType[probas.size()];
    for(unsigned int i = 0; i < probas.size(); i++)
    {
        GaussianFunctionType::MeanVectorType mu = (this->m_GaussianModel[i])->GetMean();
        (xi[i]).SetSize(mu.Size(),1);
        for(unsigned int j = 0; j < mu.Size(); j++)
        {
            (xi[i])(j,0) = mu[j];
        }
    }

    GaussianFunctionType::CovarianceMatrixType x,xT;
    for(histoIt = this->m_OriginalJointHistogram.begin(); histoIt != this->m_OriginalJointHistogram.end(); ++histoIt)
    {
        //We set the intensities of the histogram in a vector
        for(unsigned int i = 0; i < histoIt->first.size(); i++)
        {
            intensities(i,0) = static_cast<double>(histoIt->first[i]);
        }

        double concentrationValue = 0.0;

        //Calculate the exponential term and the minimum of them
        for(unsigned int i = 0; i < probas.size(); i++)
        {
            x = intensities - xi[i];
            xT = (x.GetTranspose());
            vnl_matrix<GaussianREMEstimator::NumericType> result = xT.GetVnlMatrix() * inverseCovariance[i] * x.GetVnlMatrix();
            probas[i] = result(0,0)/2;
            concentrationValue += this->m_Alphas[i] * std::exp(-probas[i]) / std::sqrt( determinantCovariance[i] );
        }

        //Add value to the GenericContainer
        // We are storing the value inside of the exponential( it will be use in expectation)
        //this->m_APosterioriProbability.insert(GenericContainer::value_type(histoIt->first,probas));
        //Fill residualmap with probabilities of the mixed gaussian (probability = constant * concentrationvalue)
        residualMap.insert(ResidualMap::value_type(std::log(concentrationValue),histoIt->first));
        numberOfPixels+=histoIt->second;
    }

    ResidualMap::iterator it = residualMap.begin();

    //number of rejected pixels
    double numberOfRejections = this->m_RejectionRatio * numberOfPixels;
    double rejected = 0;
    this->m_JointHistogram = this->m_OriginalJointHistogram;

    for(it=residualMap.begin(); it != residualMap.end(); ++it)
    {
        if(rejected >= numberOfRejections)
            break;
        histoIt = this->m_JointHistogram.find(it->second);
        if(histoIt == this->m_JointHistogram.end())
        {
            return false;
        }
        double actual = histoIt->second;
        if(actual+rejected >= numberOfRejections)
        {
            //We pass the limit...we get only some points of this Intensities
            histoIt->second = actual+rejected-numberOfRejections;
            break;
        }
        else
        {
            //We don't pass the limit... we eliminate this Intensities
            this->m_JointHistogram.erase(histoIt);
            rejected += actual;
        }
    }

    delete[] determinantCovariance;
    determinantCovariance = NULL;

    delete[] xi;
    xi = NULL;

    return true;
}

template <typename TInputImage, typename TMaskImage>
int GaussianREMEstimator<TInputImage,TMaskImage>
::PrintSolution(std::vector<double> alphas, std::vector<GaussianFunctionType::Pointer> model)
{
    unsigned int nbTissus = alphas.size();
    unsigned int nbModalities = (model[0])->GetMean().Size();
    for(unsigned int i = 0; i < nbTissus; i++)
    {
        std::cout << "* Class: " << i << std::endl;
        std::cout << "  Alpha: " << alphas[i] << std::endl;
        GaussianFunctionType::MeanVectorType mu = (model[i])->GetMean();
        std::cout << "  Mean: " << std::endl;
        for(unsigned int j = 0; j < nbModalities; j++)
        {
            std::cout << "    " << mu[j] << std::endl;
        }

        GaussianFunctionType::CovarianceMatrixType covar = (model[i])->GetCovariance();
        std::cout << "  covar: " << std::endl;
        for(unsigned int k = 0; k < nbModalities; k++)
        {
            std::cout << "    " << covar(k,0) << " " << covar(k,1) << " " << covar(k,2) << std::endl;
        }
    }
    return 0;
}

template <typename TInputImage, typename TMaskImage>
void GaussianREMEstimator<TInputImage,TMaskImage>::Update()
{
    this->createJointHistogram();

    this->m_OriginalJointHistogram = this->m_JointHistogramInitial;

    unsigned int iter = 0; //number of current iterations
    double distance = 0.0;

    this->m_Likelihood = 0.0;

    if( this->m_StremMode )
        this->m_JointHistogram = this->m_OriginalJointHistogram;
    else
    {
        if( !this->concentration() )
        {
            this->m_Likelihood = 0.0;
            return;
        }
    }

    do
    {
        int emIter = 0; // iterations between concentration steps
        do
        {
            this->m_Likelihood = this->expectation();

            if( this->m_Likelihood >= 0.0 )
            {
                this->m_Likelihood = 0.0;
                return;
            }

            std::vector<GaussianFunctionType::Pointer> newModel;
            std::vector<double> newAlphas;
            if( !this->maximization(newModel,newAlphas) )
            {
                this->m_Likelihood = 0.0;
                return;
            }

            distance = this->computeDistance(newModel);

            this->m_GaussianModel = newModel;
            this->m_Alphas = newAlphas;

            if( this->m_Verbose )
            {
                std::cout << "current iteration / max iteration between concentration step: "<< emIter << " / " << this->m_MaxIterationsConc << std::endl;
                std::cout << "current distance / minimum distance: "<< distance << " / " << this->m_ModelMinDistance << std::endl;
                std::cout << "Current model: " << std::endl;
                this->PrintSolution(this->m_Alphas, this->m_GaussianModel);
                std::cout << std::endl;

                double ratio = this->m_ModelMinDistance / distance;
                if( ratio > 1 )
                {
                    ratio = 1;
                }
                this->UpdateProgress(ratio);
                std::cout << std::endl;
            }

            emIter++;

        }while((distance > this->m_ModelMinDistance) && emIter < this->m_MaxIterationsConc);

        if( !this->concentration() )
        {
            this->m_Likelihood = 0.0;
            return;
        }

        iter++;
    }while((distance > this->m_ModelMinDistance) && iter < this->m_MaxIterations);

    this->m_Likelihood = this->expectation();
}

}
