#include "animaGaussianEMEstimator.h"

namespace anima
{

template <typename TInputImage, typename TMaskImage>
void GaussianEMEstimator<TInputImage,TMaskImage>::SetMask(const TMaskImage* mask)
{
    this->SetNthInput(0, const_cast<TMaskImage*>(mask));
}

template <typename TInputImage, typename TMaskImage>
void GaussianEMEstimator<TInputImage,TMaskImage>::SetInputImage1(const TInputImage* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInputImage*>(image));
    m_IndexImage1=m_NbInputs;
    m_NbInputs++;
}

template <typename TInputImage, typename TMaskImage>
void GaussianEMEstimator<TInputImage,TMaskImage>::SetInputImage2(const TInputImage* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInputImage*>(image));
    m_IndexImage2=m_NbInputs;
    m_NbInputs++;
}

template <typename TInputImage, typename TMaskImage>
void GaussianEMEstimator<TInputImage,TMaskImage>::SetInputImage3(const TInputImage* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInputImage*>(image));
    m_IndexImage3=m_NbInputs;
    m_NbInputs++;
}

template <typename TInputImage, typename TMaskImage>
void GaussianEMEstimator<TInputImage,TMaskImage>::SetInputImage4(const TInputImage* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInputImage*>(image));
    m_IndexImage4=m_NbInputs;
    m_NbInputs++;
}

template <typename TInputImage, typename TMaskImage>
void GaussianEMEstimator<TInputImage,TMaskImage>::SetInputImage5(const TInputImage* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInputImage*>(image));
    m_IndexImage5=m_NbInputs;
    m_NbInputs++;
}

template <typename TInputImage, typename TMaskImage>
typename TMaskImage::ConstPointer GaussianEMEstimator<TInputImage,TMaskImage>::GetMask()
{
    return static_cast< const TMaskImage * >
            ( this->itk::ProcessObject::GetInput(0) );
}

template <typename TInputImage, typename TMaskImage>
typename TInputImage::ConstPointer GaussianEMEstimator<TInputImage,TMaskImage>::GetInputImage1()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage1) );
}

template <typename TInputImage, typename TMaskImage>
typename TInputImage::ConstPointer GaussianEMEstimator<TInputImage,TMaskImage>::GetInputImage2()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage2) );
}

template <typename TInputImage, typename TMaskImage>
typename TInputImage::ConstPointer GaussianEMEstimator<TInputImage,TMaskImage>::GetInputImage3()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage3) );
}

template <typename TInputImage, typename TMaskImage>
typename TInputImage::ConstPointer GaussianEMEstimator<TInputImage,TMaskImage>::GetInputImage4()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage4) );
}

template <typename TInputImage, typename TMaskImage>
typename TInputImage::ConstPointer GaussianEMEstimator<TInputImage,TMaskImage>::GetInputImage5()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage5) );
}

template <typename TInputImage, typename TMaskImage>
void GaussianEMEstimator<TInputImage,TMaskImage>::createJointHistogram()
{
    m_ImagesVector.clear();
    m_JointHistogramInitial.clear();

    if(m_IndexImage1 < m_nbMaxImages){m_ImagesVector.push_back(this->GetInputImage1());}
    if(m_IndexImage2 < m_nbMaxImages){m_ImagesVector.push_back(this->GetInputImage2());}
    if(m_IndexImage3 < m_nbMaxImages){m_ImagesVector.push_back(this->GetInputImage3());}
    if(m_IndexImage4 < m_nbMaxImages){m_ImagesVector.push_back(this->GetInputImage4());}
    if(m_IndexImage5 < m_nbMaxImages){m_ImagesVector.push_back(this->GetInputImage5());}

    unsigned int histoDimension = m_ImagesVector.size();
    std::vector<InputConstIteratorType> ImagesVectorIt;
    for ( unsigned int i = 0; i < m_ImagesVector.size(); i++ )
    {
        InputConstIteratorType It(m_ImagesVector[i],m_ImagesVector[i]->GetLargestPossibleRegion() );
        ImagesVectorIt.push_back(It);
    }
    MaskConstIteratorType MaskIt (this->GetMask(), this->GetMask()->GetLargestPossibleRegion() );
    Histogram::iterator it;
    while (!MaskIt.IsAtEnd())
    {
        if(MaskIt.Get()!=0)
        {
            Intensities value(histoDimension);
            for(unsigned int m = 0; m < histoDimension; m++ )
            {
                if(static_cast<MeasureType>(ImagesVectorIt[m].Get()) < 0 )
                    value[m] = 0;
                if(static_cast<MeasureType>(ImagesVectorIt[m].Get()) > static_cast<double>(std::numeric_limits<MeasureType>::max()))
                    value[m]= std::numeric_limits<MeasureType>::max();
                else
                    value[m] = static_cast<MeasureType>(ImagesVectorIt[m].Get());
            }
            it = m_JointHistogramInitial.find(value);
            if(it == m_JointHistogramInitial.end())
            {
                m_JointHistogramInitial.insert(Histogram::value_type(value,1));
            }
            else
            {
                (it->second)++;
            }
        }
        for ( unsigned int i = 0; i < histoDimension; i++ )
        {
            ++ImagesVectorIt[i];
        }
        ++MaskIt;
    }
}

template <typename TInputImage, typename TMaskImage>
double GaussianEMEstimator<TInputImage,TMaskImage>::expectation()
{
    unsigned int nbClasses = m_GaussianModel.size();
    this->m_APosterioriProbability.clear();

    //1. We calculate the inverse of the covariance and the determinant;
    GaussianFunctionType::CovarianceMatrixType *inverseCovariance = new GaussianFunctionType::CovarianceMatrixType[this->m_GaussianModel.size()];
    double *determinantCovariance = new double[nbClasses];
    for(unsigned int i = 0 ; i < nbClasses; i++)
    {
        GaussianFunctionType::CovarianceMatrixType covar = (m_GaussianModel[i])->GetCovariance();

        determinantCovariance[i] = vnl_determinant(covar.GetVnlMatrix());
        if(std::abs(determinantCovariance[i]) < 1e-12)
        {
            if(inverseCovariance != NULL)
            {
                delete[] inverseCovariance;
                inverseCovariance=NULL;
            }
            if(determinantCovariance != NULL)
            {
                delete[] determinantCovariance;
                determinantCovariance= NULL;
            }
            return 1.0;
        }
        inverseCovariance[i] = covar.GetInverse();
    }

    //2. We calculate the a posterirori probability
    Histogram::iterator histoIt;
    GaussianFunctionType::CovarianceMatrixType intensities(this->m_JointHistogram.begin()->first.size(),1);
    std::vector<double> probas(nbClasses);

    GaussianFunctionType::CovarianceMatrixType *xi = new GaussianFunctionType::CovarianceMatrixType[nbClasses];
    for(unsigned int i = 0; i < nbClasses; i++)
    {
        GaussianFunctionType::MeanVectorType mu = (m_GaussianModel[i])->GetMean();
        for(unsigned int j = 0; j < mu.Size(); j++)
        {
            (xi[i]).SetSize(mu.Size(),1);
            (xi[i])(j,0) = mu[j];
        }
    }

    GaussianFunctionType::CovarianceMatrixType x,xT,result;

    for(histoIt = this->m_JointHistogram.begin(); histoIt != this->m_JointHistogram.end(); ++histoIt)
    {
        //We set the intensities of the histogram in a vector
        for(unsigned int i = 0; i < histoIt->first.size(); i++)
        {
            intensities(i,0) = static_cast<double>(histoIt->first[i]);
        }

        // Calculate probability a posteriori for each pixel
        // To eliminate problems with too small numbers we are going to substract in the exponetial
        // the minimum found to at least have one "significant" value (equivalent to multiply the whole for a constant)
        // Afterwards the a posteriory probability is normalize so this constant is eliminated
        double minExpoTerm = 1e10;
        double sumProba = 0.0;

        //Calculate the exponential term and the minimum of them
        for(unsigned int i = 0; i < probas.size(); i++)
        {
            x = intensities - xi[i];
            xT = x.GetTranspose();
            result = xT * inverseCovariance[i] * x;
            probas[i] = result(0,0);
            if( minExpoTerm > probas[i])
            {
                minExpoTerm = probas[i];
            }
        }

        //Calculate numerator and denominator of a posteriori probability
        for(unsigned int i = 0; i < probas.size(); i++)
        {
            // Constant * conditional probability of pixel knwoing gaussian i, last value multiply by the gaussian proportion
            probas[i] = m_Alphas[i] * ((std::exp(0.5* (minExpoTerm - probas[i])) / std::sqrt(std::fabs(determinantCovariance[i]))));
            // adding all classes to have the denominator
            sumProba += probas[i];
        }

        // This division eliminate the constant ment before
        for(unsigned int i = 0; i < probas.size(); i++)
        {
            probas[i] /= sumProba;
        }

        //Add probability to the GenericContainer
        this->m_APosterioriProbability.insert(GenericContainer::value_type(histoIt->first, probas));
    }

    double likelihoodValue = this->likelihood(inverseCovariance, determinantCovariance);

    // memory free
    if(inverseCovariance != NULL)
    {
        delete[] inverseCovariance;
        inverseCovariance = NULL;
    }
    if(determinantCovariance != NULL)
    {
        delete[] determinantCovariance;
        determinantCovariance = NULL;
    }

    delete[] xi;
    xi = NULL;

    return likelihoodValue;
}

template <typename TInputImage, typename TMaskImage>
bool GaussianEMEstimator<TInputImage,TMaskImage>::maximization(std::vector<GaussianFunctionType::Pointer>  &newModel, std::vector<double> &newAlphas)
{
    GenericContainer::iterator it;
    Histogram::iterator histoIt;

    unsigned int numberOfClasses = m_GaussianModel.size();
    unsigned int dimensions = this->m_JointHistogram.begin()->first.size();
    double numberOfPixels = 0;

    // initializations for estimations
    double *mixedProportions = new double[numberOfClasses];
    double **means = new double*[numberOfClasses];
    for (unsigned int i = 0; i < numberOfClasses; ++i)
    {
        means[i] = new double [dimensions];
    }

    std::vector<GaussianFunctionType::CovarianceMatrixType> covariances(numberOfClasses, GaussianFunctionType::CovarianceMatrixType(dimensions,dimensions));

    for(unsigned int i = 0; i < numberOfClasses;i++)
    {
        mixedProportions[i] = 0.0;
        for(unsigned int j = 0; j < dimensions; j++)
        {
            means[i][j] = 0.0;
            for(unsigned int k = 0; k < dimensions; k++)
                covariances[i](j,k)=0.0;
        }
    }

    //Mixing proportions and gaussian means
    for(it = this->m_APosterioriProbability.begin(),histoIt = this->m_JointHistogram.begin(); it != this->m_APosterioriProbability.end(); ++it, ++histoIt) //for all histogram
    {
        numberOfPixels += static_cast<double>(histoIt->second);//counting all pixels
        for(unsigned int i = 0; i < numberOfClasses; i++)
        {
            mixedProportions[i] += it->second[i] * histoIt->second; //mixedProportions [A posteriori probability] * [occurrences]

            for(unsigned int j = 0; j < dimensions; j++)
                means[i][j] += it->second[i] * histoIt->second * histoIt->first[j]; // means: [A posteriori probability] * [occurrences]*[intensity]
        }
    }

    for(unsigned int i = 0; i < numberOfClasses; i++)
    {
        // normalization of means by sum( [A posteriori probability] * [occurrences])
        for(unsigned int j = 0; j < dimensions; j++)
            means[i][j] /= mixedProportions[i];
    }

    // Covariance matrix for gaussians
    for(it = this->m_APosterioriProbability.begin(), histoIt = this->m_JointHistogram.begin() ; it != this->m_APosterioriProbability.end(); ++it, ++histoIt) //for all histogram
    {
        for(unsigned int i = 0; i < numberOfClasses; i++)
            for(unsigned int j = 0; j < dimensions; j++)
                for(unsigned int k = j; k < dimensions; k++)
                    covariances[i](j,k) += it->second[i] * histoIt->second * ( histoIt->first[j] - means[i][j] ) * ( histoIt->first[k] - means[i][k] ); //[post proba] [occurrences] ([intensity]-[mean])^2
    }

    double mix_pro;
    for(unsigned int i = 0; i < numberOfClasses; i++)
    {
        for(unsigned int j = 0; j < dimensions; j++)
        {
            mix_pro = mixedProportions[i];
            covariances[i](j,j) /= mixedProportions[i];
            for(unsigned int k = j+1; k < dimensions; k++)
            {
                covariances[i](j,k) /= mixedProportions[i];
                covariances[i](k,j) = covariances[i](j,k);
            }
        }
        mixedProportions[i] /= static_cast<double>(numberOfPixels); // normalization of proportions by [numberOfPixels]
    }


    //storing values in an appropiate class
    newModel.clear();
    int *sort = new int[numberOfClasses]; //sorting in increasing order the means[0]
    for (unsigned int i = 0; i < numberOfClasses;i++)
    {
        sort[i] =-1;
        double minValue = 1e32;
        for(unsigned int j = 0; j < numberOfClasses;j++)
        {
            bool used = false;
            //we check j wasn't used yet
            for(unsigned int k = 0; k < i; k++)
            {
                if(j == static_cast<unsigned int>(sort[k]))
                {
                    used=true;
                    break;
                }
            }
            // if not used we get the min
            if(!used && means[j][0] < minValue)
            {
                minValue = means[j][0];
                sort[i] = j;
            }
        }
    }

    // storing classes in first image order
    for (unsigned int i = 0; i < numberOfClasses;i++)
    {
        if(sort[i]==-1)
        {
            return false;
        }

        GaussianFunctionType::MeanVectorType mu(dimensions);
        for(unsigned int j = 0; j < dimensions; j++)
        {
            mu[j] = means[sort[i]][j];
        }

        GaussianFunctionType::Pointer tmp = GaussianFunctionType::New();
        tmp->SetMean(mu);
        tmp->SetCovariance(covariances[sort[i]]);
        newAlphas.push_back(mixedProportions[sort[i]]);
        newModel.push_back(tmp);
    }

    delete[] sort;
    for (unsigned int i = 0; i < numberOfClasses; ++i)
    {
        delete[] (means[i]);
    }
    delete[] mixedProportions;

    return true;
}

template <typename TInputImage, typename TMaskImage>
void GaussianEMEstimator<TInputImage,TMaskImage>::Update()
{
    this->createJointHistogram();
    this->m_JointHistogram = this->m_JointHistogramInitial;
    unsigned int iter = 0; //number of current iterations
    double distance = 0.0;
    m_Likelihood = 0.0;

    do
    {
        m_Likelihood = this->expectation();
        if( m_Likelihood >= 0.0 )
        {
            m_Likelihood = 0.0;
            return;
        }

        std::vector<GaussianFunctionType::Pointer> newModel;
        std::vector<double> newAlphas;

        if( !this->maximization(newModel,newAlphas) )
        {
            m_Likelihood = 0.0;
            return;
        }
        distance = this->computeDistance(newModel);
        m_GaussianModel = newModel;
        m_Alphas = newAlphas;
        iter++;

    }while((distance > this->m_ModelMinDistance) && iter < this->m_MaxIterations);
    m_Likelihood = this->expectation();
}

template <typename TInputImage, typename TMaskImage>
double GaussianEMEstimator<TInputImage,TMaskImage>::likelihood(GaussianFunctionType::CovarianceMatrixType *invCovariance, double *detCovariance)
{
    double likelihoodValue = 0.0;
    unsigned int nbClasses = m_GaussianModel.size();
    GaussianFunctionType::CovarianceMatrixType *inverseCovariance = NULL;
    double *determinantCovariance = NULL;
    if(invCovariance != NULL && detCovariance !=NULL)
    {
        inverseCovariance = invCovariance;
        determinantCovariance = detCovariance;
    }
    else
    {
        inverseCovariance = new GaussianFunctionType::CovarianceMatrixType[nbClasses];
        determinantCovariance = new double[nbClasses];

        //1. We calculate covariance inverse and determinant
        for(unsigned int i = 0 ; i < nbClasses; i++)
        {
            GaussianFunctionType::CovarianceMatrixType covar = (m_GaussianModel[i])->GetCovariance();
            determinantCovariance[i] = vnl_determinant(covar.GetVnlMatrix());
            if(std::fabs(determinantCovariance[i]) < 1e-9)
            {
                return likelihoodValue;
            }
            inverseCovariance[i] = covar.GetInverse();
        }
    }

    //2. We calculate the a posteriori probability
    Histogram::iterator histoIt;
    GenericContainer::iterator it = this->m_APosterioriProbability.begin();

    GaussianFunctionType::CovarianceMatrixType intensities(this->m_JointHistogram.begin()->first.size(),1);

    for(histoIt = this->m_JointHistogram.begin(); histoIt != this->m_JointHistogram.end(); ++histoIt, ++it)
    {
        //We set the intensities of the histogram in a vector
        for(unsigned int i = 0; i < histoIt->first.size(); i++)
        {
            intensities(i,0) = static_cast<double> (histoIt->first[i]);
        }

        unsigned int maxIndex = 0;
        double maxPostProba = 0.0;
        //we look for the max post proba to resolve de ecuation
        for(unsigned int i = 0; i < nbClasses; i++)
        {
            if(it->second[i] > maxPostProba)
            {
                maxPostProba = it->second[i];
                maxIndex = i;
            }
        }

        GaussianFunctionType::CovarianceMatrixType x,xT,result;
        GaussianFunctionType::MeanVectorType mu = (m_GaussianModel[maxIndex])->GetMean();
        for(unsigned int i = 0; i < mu.Size(); i++)
        {
            x.SetSize(mu.Size(),1);
            x(i,0) = mu[i];
        }
        x = intensities - x;
        xT = (x.GetTranspose());
        result = xT * inverseCovariance[maxIndex] * x;
        double proba = result(0,0);

        likelihoodValue += histoIt->second * ( -proba/2.0 - std::log(std::sqrt(pow(2*M_PI,static_cast<int>(this->m_JointHistogram.begin()->first.size())) *
                                                                              std::fabs(determinantCovariance[maxIndex])))
                                              + std::log(m_Alphas[maxIndex]/maxPostProba));
    }

    //Cleaning pointers if created in this function
    if((invCovariance == NULL)&&(inverseCovariance!=NULL))
    {
        delete[] inverseCovariance;
        inverseCovariance=NULL;
    }
    if((detCovariance ==NULL)&&(determinantCovariance!=NULL))
    {
        delete[] determinantCovariance;
        determinantCovariance=NULL;
    }

    return likelihoodValue;
}


template <typename TInputImage, typename TMaskImage>
double GaussianEMEstimator<TInputImage,TMaskImage>::computeDistance(std::vector<GaussianFunctionType::Pointer> &newModel)
{
    unsigned int nbClasses = m_GaussianModel.size();

    if(newModel.size() != nbClasses)
        return 1e9; //arbitrary high value

    double criterium = 0.0;

    unsigned int numberOfGaussians = (newModel[0])->GetMean().Size();

    for(unsigned int i = 0; i < m_GaussianModel.size(); i++)
    {
        GaussianFunctionType::MeanVectorType newMu = (m_GaussianModel[i])->GetMean();
        GaussianFunctionType::MeanVectorType oldMu = (newModel[i])->GetMean();

        if(newMu.Size()!=oldMu.Size())
            return 1e9;

        for(unsigned int j = 0; j < numberOfGaussians; j++)
            criterium += std::fabs( oldMu[j] - newMu[j] );
    }

    return criterium/(static_cast<double>(nbClasses * numberOfGaussians));
}

}
