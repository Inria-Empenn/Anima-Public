#include "animaRandomInitializer.h"

namespace anima
{

void RandomInitializer::Update()
{
    if(this->minValues.size() != this->maxValues.size())
    {
        std::cerr << " --error: minValues and maxValues don't have the same size" <<std::endl;
        return;
    }
    if(this->minValues.size()==0)
    {
        std::cerr << " --error: Initialization vector was empty" <<std::endl;
        return;
    }
    if(m_NbGaussian==0)
    {
        std::cerr << " --error: number of gaussian is null" <<std::endl;
        return;
    }

    m_Alphas.clear();
    m_GaussianModel.clear();
    m_DimensionGaussian = this->maxValues.size();

    //Alpha's initialization
    double *alphas = new double[m_NbGaussian];
    double norm = 0.0;
    unsigned int i = 0;

    while( i < m_NbGaussian)
    {
        double f = (static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX));
        if(f!=0)
        {
            alphas[i] = f;
            norm += f;
            i++;
        }
    }

    for( i = 0; i < m_NbGaussian; i++)
    {
        m_Alphas.push_back(alphas[i]/norm);
        GaussianFunctionType::Pointer distribution = this->randomDistribution();
        m_GaussianModel.push_back(distribution);
    }

    delete[] alphas;
    alphas = NULL;
}

itk::Statistics::GaussianMembershipFunction<itk::VariableLengthVector<double> >::Pointer RandomInitializer::randomDistribution()
{
    GaussianFunctionType::Pointer distribution = GaussianFunctionType::New();
    GaussianFunctionType::MeanVectorType mean( m_DimensionGaussian );
    mean.Fill( 0.0 );
    GaussianFunctionType::CovarianceMatrixType cov;
    cov.SetSize( m_DimensionGaussian, m_DimensionGaussian );
    cov.Fill( 0.0 );

    for(unsigned int j = 0; j < m_DimensionGaussian; j++)
    {
        mean[j] = randUniform( this->minValues[j],this->maxValues[j]);
        cov[j][j] = randUniform(this->minValues[j],this->maxValues[j])+0.001; //we don't want a 0.0 variance
    }

    distribution->SetMean( mean );
    distribution->SetCovariance( cov );

    return distribution;
}

double RandomInitializer::randUniform(double min,double max)
{
    double value = ( static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX) ); // uniform variable between 0.0 et 1.000
    return value*(max-min) + min;
}


}
