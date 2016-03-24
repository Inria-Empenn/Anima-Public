#include "animaAtlasInitializer.h"

namespace anima
{

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
void AtlasInitializer<TInputImage, TMaskImage, TAtlasImage>::SetMask(const TMaskImage* mask)
{
    this->SetNthInput(0, const_cast<TMaskImage*>(mask));
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
void AtlasInitializer<TInputImage, TMaskImage, TAtlasImage>::SetInputImage1(const TInputImage* image)
{
    this->SetNthInput(1, const_cast<TInputImage*>(image));
    m_ImagesVector.push_back(image);
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
void AtlasInitializer<TInputImage, TMaskImage, TAtlasImage>::SetInputImage2(const TInputImage* image)
{
    this->SetNthInput(2, const_cast<TInputImage*>(image));
    m_ImagesVector.push_back(image);
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
void AtlasInitializer<TInputImage, TMaskImage, TAtlasImage>::SetInputImage3(const TInputImage* image)
{
    this->SetNthInput(3, const_cast<TInputImage*>(image));
    m_ImagesVector.push_back(image);
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
void AtlasInitializer<TInputImage, TMaskImage, TAtlasImage>::SetAtlasImage1(const TAtlasImage* image)
{
    this->SetNthInput(4, const_cast<TAtlasImage*>(image));
    m_AtlasVector.push_back(image);
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
void AtlasInitializer<TInputImage, TMaskImage, TAtlasImage>::SetAtlasImage2(const TAtlasImage* image)
{
    this->SetNthInput(5, const_cast<TAtlasImage*>(image));
    m_AtlasVector.push_back(image);
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
void AtlasInitializer<TInputImage, TMaskImage, TAtlasImage>::SetAtlasImage3(const TAtlasImage* image)
{
    this->SetNthInput(6, const_cast<TAtlasImage*>(image));
    m_AtlasVector.push_back(image);
}


template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
typename TMaskImage::ConstPointer AtlasInitializer<TInputImage, TMaskImage, TAtlasImage>::GetMask()
{
    return static_cast< const TMaskImage * >
            ( this->itk::ProcessObject::GetInput(0) );
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
typename TInputImage::ConstPointer AtlasInitializer<TInputImage, TMaskImage, TAtlasImage>::GetInputImage1()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(1) );
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
typename TInputImage::ConstPointer AtlasInitializer<TInputImage, TMaskImage, TAtlasImage>::GetInputImage2()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(2) );
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
typename TInputImage::ConstPointer AtlasInitializer<TInputImage, TMaskImage, TAtlasImage>::GetInputImage3()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(3) );
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
typename TAtlasImage::ConstPointer AtlasInitializer<TInputImage, TMaskImage, TAtlasImage>::GetAtlasImage1()
{
    return static_cast< const TAtlasImage * >
            ( this->itk::ProcessObject::GetInput(4) );
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
typename TAtlasImage::ConstPointer AtlasInitializer<TInputImage, TMaskImage, TAtlasImage>::GetAtlasImage2()
{
    return static_cast< const TAtlasImage * >
            ( this->itk::ProcessObject::GetInput(5) );
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
typename TAtlasImage::ConstPointer AtlasInitializer<TInputImage, TMaskImage, TAtlasImage>::GetAtlasImage3()
{
    return static_cast< const TAtlasImage * >
            ( this->itk::ProcessObject::GetInput(6) );
}



template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
void AtlasInitializer<TInputImage, TMaskImage, TAtlasImage>::Update()
{
    const unsigned int numberOfClasses = m_ImagesVector.size();

    /*--------------------Creates iterators --------------------*/
    MaskConstIteratorType MaskIt(this->GetMask(),this->GetMask()->GetLargestPossibleRegion());

    std::vector<InputConstIteratorType> ImagesVectorVectorIt;
    for ( unsigned int i = 0; i < numberOfClasses; i++ )
    {
        InputConstIteratorType It(m_ImagesVector[i],m_ImagesVector[i]->GetLargestPossibleRegion() );
        ImagesVectorVectorIt.push_back(It);
    }

    std::vector<AtlasConstIteratorType> AtlasVectorVectorIt;
    for ( unsigned int i = 0; i < numberOfClasses; i++ )
    {
        AtlasConstIteratorType It(m_AtlasVector[i],m_AtlasVector[i]->GetLargestPossibleRegion() );
        AtlasVectorVectorIt.push_back(It);
    }

    std::vector<GaussianFunctionType::MeanVectorType> means(numberOfClasses,GaussianFunctionType::MeanVectorType( numberOfClasses ));
    std::vector<DoubleVariableSizeMatrixType> covariances(numberOfClasses,DoubleVariableSizeMatrixType(numberOfClasses,numberOfClasses));

    m_Alphas.resize(numberOfClasses,0);

    /*--------------------Initialisation alpha, mean and covariances --------------------*/
    for(unsigned int c = 0; c < numberOfClasses; c++)
    {
        for(unsigned int i = 0; i < numberOfClasses; i++)
        {
            means[c][i] = 0.0;
            for(unsigned int j = 0; j < numberOfClasses; j++)
                covariances[c](i,j) = 0.0;
        }
    }

    /*--------------------Computes Means --------------------*/
    MaskIt.GoToBegin();
    for ( unsigned int k = 0; k < numberOfClasses; k++ )
    {
        AtlasVectorVectorIt[k].GoToBegin();
        ImagesVectorVectorIt[k].GoToBegin();
    }

    while(!MaskIt.IsAtEnd())
    {
        if ( MaskIt.Get() != 0 )
        {
            for(unsigned int c = 0; c < numberOfClasses; c++)
            {
                //Getting factors of normalization for apriori m_ImagesVector
                m_Alphas[c]+=static_cast<double>(AtlasVectorVectorIt[c].Get());

                //Mean values
                for(unsigned int m = 0; m < numberOfClasses; m++)
                {
                    (means[c])[m] += AtlasVectorVectorIt[c].Get() * ImagesVectorVectorIt[m].Get();
                }
            }
        }
        ++MaskIt;
        for ( unsigned int k = 0; k < numberOfClasses; k++ )
        {
            ++AtlasVectorVectorIt[k];
            ++ImagesVectorVectorIt[k];
        }
    }

    //normalization of mean
    for(unsigned int c = 0; c < numberOfClasses; c++)
    {
        for(unsigned int m = 0; m < numberOfClasses; m++)
        {
            (means[c])[m]/=m_Alphas[c];
        }
    }


    /*-------------------- Compute Covariances --------------------*/
    MaskIt.GoToBegin();
    for ( unsigned int i = 0; i < numberOfClasses; i++ )
    {
        ImagesVectorVectorIt[i].GoToBegin();
        AtlasVectorVectorIt[i].GoToBegin();
    }

    //Covariance
    DoubleVariableSizeMatrixType x(numberOfClasses,1);
    DoubleVariableSizeMatrixType xT (1,numberOfClasses);
    DoubleVariableSizeMatrixType meanT(1,numberOfClasses);
    DoubleVariableSizeMatrixType meanTTrans(numberOfClasses,1);

    while(!MaskIt.IsAtEnd())
    {
        if ( MaskIt.Get() != 0 )
        {
            for(unsigned int c = 0; c < numberOfClasses; c++)
            {
                for(unsigned int m = 0; m < numberOfClasses; m++)
                {
                    x(m,0) = ImagesVectorVectorIt[m].Get();
                    xT(0,m) = x(m,0);
                    meanT(0,m) = (means[c])[m];
                    meanTTrans(m,0) = (means[c])[m];
                }
                covariances[c] += ( (x-meanTTrans) * (xT-meanT) ) * (static_cast<double>(AtlasVectorVectorIt[c].Get()) );
            }
        }
        ++MaskIt;
        for ( unsigned int i = 0; i < numberOfClasses; i++ )
        {
            ++ImagesVectorVectorIt[i];
            ++AtlasVectorVectorIt[i];
        }
    }


    /*-------------------- Compute Covariances Normalisation --------------------*/
    for(unsigned int c = 0; c < numberOfClasses; c++)
    {
        //Covariance normalization
        for(unsigned int m = 0; m < numberOfClasses; m++)
        {
            for(unsigned int n = 0; n < numberOfClasses; n++)
            {
                covariances[c](m,n) /= m_Alphas[c];
            }
        }
    }

    double total=0.0;
    for(unsigned int c = 0; c < numberOfClasses; c++)
    {
        total+=m_Alphas[c];
    }

    for(unsigned int c = 0; c < numberOfClasses; c++)
    {
        GaussianFunctionType::Pointer gaussian = GaussianFunctionType::New();
        gaussian->SetMean(means[c]);
        gaussian->SetCovariance(covariances[c]);
        m_GaussianModel.push_back(gaussian);
        m_Alphas[c]/=total; // alpha normalization
    }

}

}
