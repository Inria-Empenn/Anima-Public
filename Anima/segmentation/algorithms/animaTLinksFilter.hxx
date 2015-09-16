#pragma once

#include "animaTLinksFilter.h"


namespace anima
{

template <typename TInput, typename TOutput>
void TLinksFilter<TInput, TOutput>::SetInputImage1(const TInput* image)
{
    this->SetNthInput(0, const_cast<TInput*>(image));
    m_NbModalities++;
    InConstIteratorType inputIterator (image, image->GetLargestPossibleRegion());
    m_imagesVectorIt.push_back(inputIterator);
    m_imagesVector.push_back(image);
}

template <typename TInput, typename TOutput>
void TLinksFilter<TInput, TOutput>::SetInputSeedSourcesMask(const TSeedMask* mask)
{
    this->SetNthInput(m_NbInputs, const_cast<TSeedMask*>(mask));
    m_IndexSourcesMask = m_NbInputs;
    m_NbInputs++;
}

template <typename TInput, typename TOutput>
void TLinksFilter<TInput, TOutput>::SetInputSeedSinksMask(const TSeedMask* mask)
{
    this->SetNthInput(m_NbInputs, const_cast<TSeedMask*>(mask));
    m_IndexSinksMask = m_NbInputs;
    m_NbInputs++;
}

template <typename TInput, typename TOutput>
void TLinksFilter<TInput, TOutput>::SetInputSeedSourcesProba(const TSeedProba* mask)
{
    this->SetNthInput(m_NbInputs, const_cast<TSeedProba*>(mask));
    m_IndexSourcesProba = m_NbInputs;
    m_NbInputs++;
}

template <typename TInput, typename TOutput>
void TLinksFilter<TInput, TOutput>::SetInputSeedSinksProba(const TSeedProba* mask)
{
    this->SetNthInput(m_NbInputs, const_cast<TSeedProba*>(mask));
    m_IndexSinksProba = m_NbInputs;
    m_NbInputs++;
}

template <typename TInput, typename TOutput>
void TLinksFilter<TInput, TOutput>::SetInputImage2(const TInput* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInput*>(image));
    m_IndexImage2 = m_NbInputs;
    m_NbInputs++;
    m_NbModalities++;
    InConstIteratorType inputIterator (image, image->GetLargestPossibleRegion());
    m_imagesVectorIt.push_back(inputIterator);
    m_imagesVector.push_back(image);
}

template <typename TInput, typename TOutput>
void TLinksFilter<TInput, TOutput>::SetInputImage3(const TInput* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInput*>(image));
    m_IndexImage3 = m_NbInputs;
    m_NbInputs++;
    m_NbModalities++;
    InConstIteratorType inputIterator (image, image->GetLargestPossibleRegion());
    m_imagesVectorIt.push_back(inputIterator);
    m_imagesVector.push_back(image);
}

template <typename TInput, typename TOutput>
void TLinksFilter<TInput, TOutput>::SetInputImage4(const TInput* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInput*>(image));
    m_IndexImage4 = m_NbInputs;
    m_NbInputs++;
    m_NbModalities++;
    InConstIteratorType inputIterator (image, image->GetLargestPossibleRegion());
    m_imagesVectorIt.push_back(inputIterator);
    m_imagesVector.push_back(image);
}

template <typename TInput, typename TOutput>
void TLinksFilter<TInput, TOutput>::SetInputImage5(const TInput* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInput*>(image));
    m_IndexImage5 = m_NbInputs;
    m_NbInputs++;
    m_NbModalities++;
    InConstIteratorType inputIterator (image, image->GetLargestPossibleRegion());
    m_imagesVectorIt.push_back(inputIterator);
    m_imagesVector.push_back(image);
}




template <typename TInput, typename TOutput>
typename TInput::ConstPointer TLinksFilter<TInput, TOutput>::GetInputImage1()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(0) );
}

template <typename TInput, typename TOutput>
itk::Image <unsigned char,3>::ConstPointer TLinksFilter<TInput, TOutput>::GetInputSeedSourcesMask()
{
    return static_cast< const TSeedMask * >
            ( this->itk::ProcessObject::GetInput(m_IndexSourcesMask) );
}

template <typename TInput, typename TOutput>
itk::Image <unsigned char,3>::ConstPointer TLinksFilter<TInput, TOutput>::GetInputSeedSinksMask()
{
    return static_cast< const TSeedMask * >
            ( this->itk::ProcessObject::GetInput(m_IndexSinksMask) );
}

template <typename TInput, typename TOutput>
itk::Image <float,3>::ConstPointer TLinksFilter<TInput, TOutput>::GetInputSeedSourcesProba()
{
    return static_cast< const TSeedProba * >
            ( this->itk::ProcessObject::GetInput(m_IndexSourcesProba) );
}

template <typename TInput, typename TOutput>
itk::Image <float,3>::ConstPointer TLinksFilter<TInput, TOutput>::GetInputSeedSinksProba()
{
    return static_cast< const TSeedProba * >
            ( this->itk::ProcessObject::GetInput(m_IndexSinksProba) );
}

template <typename TInput, typename TOutput>
typename TInput::ConstPointer TLinksFilter<TInput, TOutput>::GetInputImage2()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage2) );
}

template <typename TInput, typename TOutput>
typename TInput::ConstPointer TLinksFilter<TInput, TOutput>::GetInputImage3()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage3) );
}
template <typename TInput, typename TOutput>
typename TInput::ConstPointer TLinksFilter<TInput, TOutput>::GetInputImage4()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage4) );
}

template <typename TInput, typename TOutput>
typename TInput::ConstPointer TLinksFilter<TInput, TOutput>::GetInputImage5()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage5) );
}



template <typename TInput, typename TOutput>
itk::DataObject::Pointer TLinksFilter<TInput, TOutput>::MakeOutput(unsigned int idx)
{
    itk::DataObject::Pointer output;

    switch ( idx )
    {
    case 0:
        output = ( TOutput::New() ).GetPointer();
        break;
    case 1:
        output = ( TOutput::New() ).GetPointer();
        break;
    default:
        std::cerr << "No output " << idx << std::endl;
        output = NULL;
        break;
    }
    return output.GetPointer();
}

template <typename TInput, typename TOutput>
TOutput* TLinksFilter<TInput, TOutput>::GetOutputSources()
{
    return dynamic_cast< TOutput * >( this->itk::ProcessObject::GetOutput(0) );
}

template <typename TInput, typename TOutput>
TOutput* TLinksFilter<TInput, TOutput>::GetOutputSinks()
{
    return dynamic_cast< TOutput * >( this->itk::ProcessObject::GetOutput(1) );
}


template <typename TInput, typename TOutput>
void
TLinksFilter<TInput, TOutput>
::GenerateData()
{
    std::cout << "Computing T-Links..." << std::endl;

    typename TOutput::Pointer output1 = this->GetOutputSources();
    output1->SetRegions(this->GetInputImage1()->GetLargestPossibleRegion());
    output1->CopyInformation(this->GetInputImage1());
    output1->Allocate();
    output1->FillBuffer(0);

    typename TOutput::Pointer output2 = this->GetOutputSinks();
    output2->SetRegions(this->GetInputImage1()->GetLargestPossibleRegion());
    output2->CopyInformation(this->GetInputImage1());
    output2->Allocate();
    output2->FillBuffer(0);

    if(this->m_TLinkMode==singleGaussianTLink && (this->GetInputSeedSourcesMask().IsNull() && this->GetInputSeedSinksMask().IsNull()))
    {
        std::cerr << "-- Error in T-Links computation: Single Gaussian mode requires sources or sinks mask" << std::endl;
        exit(-1);
    }
    if(this->m_TLinkMode==stremTLink && (this->GetInputSeedSourcesProba().IsNull() && this->GetInputSeedSinksProba().IsNull()))
    {
        std::cerr << "-- Error in T-Links computation: Strem mode requires sources or sinks probabilties" << std::endl;
        exit(-1);
    }

    m_NbModalities = m_imagesVectorIt.size();

    switch(m_TLinkMode)
    {
    case singleGaussianTLink:
        computeSingleGaussian();
        break;
    case stremTLink:
        computeStrem();
        break;
    default:
        return;
    }

    OutRegionIteratorType outIterator1(output1, output1->GetLargestPossibleRegion());
    OutRegionIteratorType outIterator2(output2, output2->GetLargestPossibleRegion());

    while (!outIterator1.IsAtEnd())
    {
        float val1 = m_Alpha * outIterator1.Get();
        outIterator1.Set(val1);

        float val2 = m_Alpha * outIterator2.Get();
        outIterator2.Set(val2);

        ++outIterator1;
        ++outIterator2;
    }
}

template <typename TInput, typename TOutput>
void
TLinksFilter<TInput, TOutput>
::computeStrem()
{
    // Just copy proba
    std::cout << "Use of strem method..." << std::endl;
    OutRegionIteratorType outIterator1(this->GetOutputSources(), this->GetOutputSources()->GetLargestPossibleRegion());
    OutRegionIteratorType outIterator2(this->GetOutputSinks(), this->GetOutputSinks()->GetLargestPossibleRegion());

    SeedProbaRegionConstIteratorType sourcesIterator(this->GetInputSeedSourcesProba(), this->GetInputSeedSourcesProba()->GetLargestPossibleRegion());
    SeedProbaRegionConstIteratorType sinksIterator(this->GetInputSeedSinksProba(), this->GetInputSeedSinksProba()->GetLargestPossibleRegion());

    while(!sourcesIterator.IsAtEnd())
    {
        outIterator1.Set(static_cast<OutputPixelType>(sourcesIterator.Get()));
        outIterator2.Set(static_cast<OutputPixelType>(sinksIterator.Get()));

        ++outIterator1;
        ++outIterator2;
        ++sourcesIterator;
        ++sinksIterator;
    }
}

template <typename TInput, typename TOutput>
void
TLinksFilter<TInput, TOutput>
::computeSingleGaussian()
{
    if(this->GetInputSeedSourcesMask().IsNotNull())
    {
        std::cout << "Computing single gaussian for sources..." << std::endl;
        computeSingleGaussianSeeds(this->GetInputSeedSourcesMask(), this->GetOutputSources(), m_MultiVarSources, this->GetInputSeedSinksMask());
    }
    if(this->GetInputSeedSinksMask().IsNotNull())
    {
        std::cout << "Computing single gaussian for sinks..." << std::endl;
        computeSingleGaussianSeeds(this->GetInputSeedSinksMask(), this->GetOutputSinks(), m_MultiVarSinks, this->GetInputSeedSourcesMask());
    }
}

template <typename TInput, typename TOutput>
void
TLinksFilter<TInput, TOutput>
::computeSingleGaussianSeeds(TSeedMask::ConstPointer seedMask, OutputImagePointer output, float multiVar, TSeedMask::ConstPointer seedMaskOpp)
{
    OutRegionIteratorType outIterator(output, output->GetLargestPossibleRegion());
    SeedMaskRegionConstIteratorType seedIterator(seedMask, seedMask->GetLargestPossibleRegion());

    double epsilon=1e-37;

    // Compute the mean and variance of each seeds for each modality
    std::vector<double> moy, var;
    moy.resize(m_NbModalities,0);
    var.resize(m_NbModalities,0);

    double nbSeeds = 0.0;
    std::vector<double> som;
    som.resize(m_NbModalities,0);

    std::vector<double> som2;
    som2.resize(m_NbModalities,0);

    for ( unsigned int k = 0; k < m_NbModalities; k++ )
    {
        m_imagesVectorIt[k].GoToBegin();
    }
    while(!seedIterator.IsAtEnd())
    {
        if(seedIterator.Get()!=0)
        {
            nbSeeds++;
            for (unsigned int m = 0; m < m_NbModalities; m++)
            {
                som [m] += m_imagesVectorIt[m].Get();
                som2[m] += std::pow(static_cast<double>(m_imagesVectorIt[m].Get()),2.0);
            }
        }
        ++seedIterator;
        for ( unsigned int k = 0; k < m_NbModalities; k++ )
        {
            ++m_imagesVectorIt[k];
        }
    }

    for (unsigned int m = 0; m < m_NbModalities; m++)
    {
        moy[m] = som[m]/nbSeeds;
        var[m] = ((som2[m]/nbSeeds) - (moy[m]*moy[m]));
        var[m]*=multiVar;
        if (var[m]==0)
        {
            var[m]=epsilon;
        }

        if(m_Verbose)
        {
            std::cout << "mean of modality " << m << ": " << moy[m] << std::endl;
            std::cout << "variance of modality " << m << ": " << var[m] << std::endl;
            std::cout << std::endl;
        }
    }

    som.clear();

    // Compute the covariances between modalities
    std::vector<double> covar;
    if (m_NbModalities > 1)
    {
        int nb_paires = 0;
        for (unsigned int m = 0; m < m_NbModalities-1; m++)
        {
            for (unsigned int n = m+1; n < m_NbModalities; n++)
            {
                nb_paires++;
            }
        }
        covar.resize(nb_paires,0);
        som.resize(nb_paires,0);

        seedIterator.GoToBegin();
        for ( unsigned int i = 0; i < m_NbModalities; i++ )
        {
            m_imagesVectorIt[i].GoToBegin();
        }

        while(!seedIterator.IsAtEnd())
        {
            if(seedIterator.Get()!=0)
            {
                int pos=0;
                for (unsigned int m = 0; m < m_NbModalities-1; m++)
                {
                    for (unsigned int n = m+1; n < m_NbModalities; n++)
                    {
                        som[pos] += (m_imagesVectorIt[m].Get()-moy[m]) * (m_imagesVectorIt[n].Get()-moy[n]);
                        pos++;
                    }
                }
            }
            ++seedIterator;
            for ( unsigned int k = 0; k < m_NbModalities; k++ )
            {
                ++m_imagesVectorIt[k];

            }
        }


        int pos = 0;
        for (unsigned int m = 0; m < m_NbModalities-1; m++)
        {
            for (unsigned int n = m+1; n < m_NbModalities; n++)
            {
                covar[pos] = som[pos]/(nbSeeds-1);
                if (covar[pos]==0)
                    covar[pos]=epsilon;
                pos++;
            }
        }
    }


    DoubleVariableSizeMatrixType covarMatrix;
    covarMatrix.SetSize(m_NbModalities,m_NbModalities);
    for (unsigned int m = 0; m < m_NbModalities; m++)
    {
        covarMatrix(m,m) = var[m];
    }
    int pos = 0;
    for (unsigned int m = 0; m < m_NbModalities-1; m++)
    {
        for (unsigned int n = m+1; n < m_NbModalities; n++)
        {
            covarMatrix(m,n) = covarMatrix(n,m) = covar[pos];
            pos++;
        }
    }

    //invert covariance matrix
    covarMatrix=covarMatrix.GetInverse();

    // Compute the proba maps
    MatrixTypeRes prob;
    DoubleVariableSizeMatrixType d;
    d.SetSize(1,m_NbModalities);

    seedIterator.GoToBegin();
    for ( unsigned int i = 0; i < m_NbModalities; i++ )
    {
        m_imagesVectorIt[i].GoToBegin();
    }

    TSeedMask::Pointer seedMaskOpp2 = TSeedMask::New();
    seedMaskOpp2->SetRegions(seedMask->GetLargestPossibleRegion());
    seedMaskOpp2->CopyInformation(seedMask);
    seedMaskOpp2->Allocate();

    SeedMaskRegionIteratorType seedOppIterator2(seedMaskOpp2, seedMaskOpp2->GetLargestPossibleRegion());
    if(seedMaskOpp.IsNull())
    {
        seedMaskOpp2->FillBuffer(0);
    }
    else
    {
        SeedMaskRegionConstIteratorType seedOppIterator(seedMaskOpp, seedMaskOpp->GetLargestPossibleRegion());
        while(!seedOppIterator.IsAtEnd())
        {
            seedOppIterator2.Set(seedOppIterator.Get());
            ++seedOppIterator2;
            ++seedOppIterator;
        }
    }

    seedIterator.GoToBegin();
    seedOppIterator2.GoToBegin();
    while(!seedIterator.IsAtEnd())
    {
        for (unsigned int m = 0; m < m_NbModalities; m++)
        {
            d(0,m) = m_imagesVectorIt[m].Get() - moy[m];
        }

        prob = d.GetVnlMatrix()*covarMatrix.GetVnlMatrix()*d.GetTranspose();

        outIterator.Set( std::exp(-0.5*prob(0,0)) );

        if(seedIterator.Get()!=0)
        {
            outIterator.Set(1.0);
        }
        if(seedOppIterator2.Get()!=0)
        {
            outIterator.Set(0.0);
        }
        ++seedOppIterator2;
        ++outIterator;
        ++seedIterator;
        for ( unsigned int k = 0; k < m_NbModalities; k++ )
        {
            ++m_imagesVectorIt[k];
        }
    }

}


} //end of namespace anima
