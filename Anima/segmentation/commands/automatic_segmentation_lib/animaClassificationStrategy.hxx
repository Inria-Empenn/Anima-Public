#include "animaClassificationStrategy.h"


namespace anima
{
template <typename TInputImage, typename TMaskImage>
void ClassificationStrategy<TInputImage,TMaskImage>::Update()
{
    if( this->m_RandomInitializer.IsNull() )
    {
        std::cerr <<" -- Error in ClassificationStrategy: RandomInitializer is null"<<std::endl;
        return;
    }

    m_ListGaussianModels.clear();
    m_ListAlphas.clear();

    m_ListGaussianModels.resize(this->m_NumberOfEstimators.size());
    m_ListAlphas.resize(this->m_NumberOfEstimators.size());

    // First iteration with RandomInitialization
    unsigned int errors = 0;

    for(unsigned int iter = 0; iter < this->m_NumberOfEstimators[0]; iter++)
    {
        this->m_RandomInitializer->Update();

        std::vector<GaussianFunctionType::Pointer> model = this->m_RandomInitializer->GetInitialization();
        std::vector<double> alphas = this->m_RandomInitializer->GetAlphas();

        this->m_Estimator->SetInitialGaussianModel(model);
        this->m_Estimator->SetInitialAlphas(alphas);
        this->m_Estimator->Update();
        double likelihood = m_Estimator->GetLikelihood();

        if( (!m_EM_Mode && likelihood > 0.0) || (m_EM_Mode && likelihood < 0.0) )
        {
            std::vector<GaussianFunctionType::Pointer> result = this->m_Estimator->GetGaussianModel();
            std::vector<double> resultAlpha = this->m_Estimator->GetAlphas();

            if(m_EM_Mode)
            {
                likelihood=-likelihood; //we change the sign to have the best solution in the first position
            }

            std::map< double, std::vector<GaussianFunctionType::Pointer> >::iterator mapIt = this->m_ListGaussianModels[0].begin();
            bool foundSame = false;

            for(mapIt = this->m_ListGaussianModels[0].begin(); mapIt != this->m_ListGaussianModels[0].end(); ++mapIt)
            {
                if(sameModel(mapIt->second, result))
                {
                    foundSame = true;
                    break;
                }

            }

            if(!foundSame)
            {
                this->m_ListGaussianModels[0].insert(std::map<double, std::vector<GaussianFunctionType::Pointer> >::value_type(likelihood,result));
                this->m_ListAlphas[0].insert(std::map<double,std::vector<double> >::value_type(likelihood,resultAlpha));
            }

        }
        else
        {
            iter--;
            errors++;
            if(errors > 10 * this->m_NumberOfEstimators[0])
                break;
        }


        double ratio = static_cast<double>(iter+1) / static_cast<double>(this->m_NumberOfEstimators[0]);
        this->UpdateProgress(ratio);
    }


    // if several estimators
    for(unsigned int step = 1; step < this->m_NumberOfEstimators.size(); step++)
    {
        if(this->m_ListGaussianModels[step-1].size() == 0)
        {
            std::cout<<"-- ERROR in ClassificationStrategy: No solution found in step " << step-1 <<std::endl;
            return;
        }

        std::map<double,std::vector<GaussianFunctionType::Pointer> >::iterator mapIt = this->m_ListGaussianModels[step-1].begin();
        for(unsigned int iter = 0 ; mapIt != this->m_ListGaussianModels[step-1].end() && iter< this->m_NumberOfEstimators[step]; iter++, ++mapIt)
        {
            this->m_Estimator->SetInitialGaussianModel(mapIt->second);

            m_Estimator->Update();
            double likelihood = m_Estimator->GetLikelihood();

            if((m_EM_Mode && likelihood < 0.0) || (!m_EM_Mode && likelihood > 0.0))
            {
                std::vector<GaussianFunctionType::Pointer> result = m_Estimator->GetGaussianModel();
                std::vector<double> resultAlpha = m_Estimator->GetAlphas();

                mapIt = this->m_ListGaussianModels[step].begin();
                bool foundSame = false;
                for(mapIt = this->m_ListGaussianModels[step].begin(); mapIt != this->m_ListGaussianModels[step].end(); ++mapIt)
                {
                    if(sameModel(mapIt->second, result))
                    {
                        foundSame=true;
                        break;
                    }
                }

                if(m_EM_Mode)
                    likelihood=-likelihood; //we change sign to have best m_Solutions first

                if(!foundSame)
                    this->m_ListGaussianModels[step].insert(std::map<double,std::vector<GaussianFunctionType::Pointer> >::value_type(likelihood,result));
                    this->m_ListAlphas[step].insert(std::map<double,std::vector<double> >::value_type(likelihood,resultAlpha));
            }
            else
            {
                //error
                iter--;
                errors++;
                if(errors > 10 * this->m_NumberOfEstimators[step])
                    break;
            }
        }
    }
}

template <typename TInputImage, typename TMaskImage>
bool ClassificationStrategy<TInputImage,TMaskImage>::GetSolutionMap(std::map< double, std::vector<GaussianFunctionType::Pointer> > &solution,
                                                                    std::map< double, std::vector<double> > &solutionAlpha, int step)
{
    if(this->m_ListGaussianModels.size() == 0 || static_cast<int>(this->m_ListGaussianModels.size()) <= step)
        return false;

    if(step < 0 )
    {
        step=this->m_ListGaussianModels.size()-1;
    }

    if(this->m_ListGaussianModels[step].size() == 0 )
        return false;

    solution = this->m_ListGaussianModels[step];
    solutionAlpha = this->m_ListAlphas[step];
    return true;
}


template <typename TInputImage, typename TMaskImage>
void ClassificationStrategy<TInputImage,TMaskImage>::SetStrategy(  std::vector< unsigned int >& ems,std::vector<unsigned int> &iters )
{
    this->m_NumberOfEstimators=ems;
    this->m_NumberOfIterations=iters;

    if(this->m_NumberOfEstimators.size()!= this->m_NumberOfIterations.size())
    {
        while( this->m_NumberOfEstimators.size() > this->m_NumberOfIterations.size() )
        {
            this->m_NumberOfEstimators.pop_back();
        }
        while( this->m_NumberOfEstimators.size() < this->m_NumberOfIterations.size() )
        {
            this->m_NumberOfIterations.pop_back();
        }
    }

}

template <typename TInputImage, typename TMaskImage>
bool ClassificationStrategy<TInputImage,TMaskImage>::sameModel( std::vector<GaussianFunctionType::Pointer> &mod1, std::vector<GaussianFunctionType::Pointer> &mod2)
{
    unsigned int NbClasses = mod1.size();
    unsigned int NbDimension = (mod1[0])->GetMean().Size();

    std::vector<int> comparisonVector(NbClasses,-1);
    unsigned int params = NbDimension + NbDimension*NbDimension;

    for(unsigned int i = 0; i < NbClasses; i++)
    {
        for(unsigned int j = 0; j < NbClasses; j++)
        {
            if((mod1[i])->GetMean().Size() != (mod2[j])->GetMean().Size())
            {
                comparisonVector[i] = -1;
                continue;
            }

            float *p1 = new float[params];
            float *p2 = new float[params];

            unsigned int t = 0;
            for(unsigned int l = 0; l < NbDimension; l++)
            {
                p1[t] = static_cast<float>((mod1[i])->GetMean()[l]);
                p2[t] = static_cast<float>((mod2[j])->GetMean()[l]);
                t++;
            }

            for(unsigned int l = 0; l < NbDimension; l++)
            {
                for(unsigned int k = l; k < NbDimension; k++)
                {
                    p1[t] = static_cast<float>((mod1[i])->GetCovariance()[l][k]);
                    p2[t] = static_cast<float>((mod2[j])->GetCovariance()[l][k]);
                    t++;
                }
            }

            double distance = 0.0;
            for(unsigned int k = 0; k < params; k++)
            {
                distance+=p1[k]-p2[k];
            }

            distance/=static_cast<double>(params);
            if(distance < 1e-2)
            {
                comparisonVector[i]=j;
                break;
            }
            delete[] p1;
            delete[] p2;
        }
    }

   for (unsigned int i=0; i < NbClasses; i++)
    {
        if(comparisonVector[i]==-1)
            return false;
        for(unsigned int j=i+1;j < NbClasses; j++)
            if( comparisonVector[i]==comparisonVector[j])
                return false;
    }

    return true;
}


}
