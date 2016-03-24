#pragma once

#include "animaImageClassifierFilter.h"

namespace anima
{

template <typename TInput, typename TMask, typename TOutput>
void ImageClassifierFilter<TInput, TMask, TOutput>::SetInputImage1(const TInput* image)
{
    this->SetNthInput(0, const_cast<TInput*>(image));
    m_ImagesVector.push_back(image);
}

template <typename TInput, typename TMask, typename TOutput>
void ImageClassifierFilter<TInput, TMask, TOutput>::SetMask(const TMask* mask)
{
    this->SetNthInput(1, const_cast<TMask*>(mask));
}

template <typename TInput, typename TMask, typename TOutput>
void ImageClassifierFilter<TInput, TMask, TOutput>::SetInputImage2(const TInput* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInput*>(image));
    m_IndexImage2 = m_NbInputs;
    m_ImagesVector.push_back(image);
    m_NbInputs++;
}

template <typename TInput, typename TMask, typename TOutput>
void ImageClassifierFilter<TInput, TMask, TOutput>::SetInputImage3(const TInput* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInput*>(image));
    m_IndexImage3 = m_NbInputs;
    m_ImagesVector.push_back(image);
    m_NbInputs++;
}
template <typename TInput, typename TMask, typename TOutput>
void ImageClassifierFilter<TInput, TMask, TOutput>::SetInputImage4(const TInput* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInput*>(image));
    m_IndexImage4 = m_NbInputs;
    m_ImagesVector.push_back(image);
    m_NbInputs++;
}

template <typename TInput, typename TMask, typename TOutput>
void ImageClassifierFilter<TInput, TMask, TOutput>::SetInputImage5(const TInput* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInput*>(image));
    m_IndexImage5 = m_NbInputs;
    m_ImagesVector.push_back(image);
    m_NbInputs++;
}



template <typename TInput, typename TMask, typename TOutput>
typename TInput::ConstPointer ImageClassifierFilter<TInput, TMask, TOutput>::GetInputImage1()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(0) );
}

template <typename TInput, typename TMask, typename TOutput>
typename TMask::ConstPointer ImageClassifierFilter<TInput, TMask, TOutput>::GetMask()
{
    return static_cast< const TMask * >
            ( this->itk::ProcessObject::GetInput(1) );
}

template <typename TInput, typename TMask, typename TOutput>
typename TInput::ConstPointer ImageClassifierFilter<TInput, TMask, TOutput>::GetInputImage2()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage2) );
}

template <typename TInput, typename TMask, typename TOutput>
typename TInput::ConstPointer ImageClassifierFilter<TInput, TMask, TOutput>::GetInputImage3()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage3) );
}

template <typename TInput, typename TMask, typename TOutput>
typename TInput::ConstPointer ImageClassifierFilter<TInput, TMask, TOutput>::GetInputImage4()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage4) );
}

template <typename TInput, typename TMask, typename TOutput>
typename TInput::ConstPointer ImageClassifierFilter<TInput, TMask, TOutput>::GetInputImage5()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage5) );
}

template <typename TInput, typename TMask, typename TOutput>
void
ImageClassifierFilter <TInput, TMask, TOutput>::WriteOutputs()
{
    if( m_OutputFilename != "" )
    {
        std::cout << "Writing classification image to: " << m_OutputFilename << std::endl;
        anima::writeImage<TOutput>(m_OutputFilename, this->GetOutput());
    }
}

template <typename TInput, typename TMask, typename TOutput>
void
ImageClassifierFilter<TInput, TMask, TOutput>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    OutputIteratorType classificationIt(this->GetOutput(), outputRegionForThread);
    MaskConstIteratorType maskIt(this->GetMask(), outputRegionForThread);

    std::vector<InputConstIteratorType > inputVectorIt;
    for ( unsigned int i = 0; i < m_ImagesVector.size(); i++ )
    {
        InputConstIteratorType It(m_ImagesVector[i], outputRegionForThread );
        inputVectorIt.push_back(It);
    }

    while(!classificationIt.IsAtEnd())
    {
        classificationIt.Set(0);
        if ( maskIt.Get() != 0 )
        {
            // fill a vector with images information
            DoubleVariableSizeMatrixType intensities(inputVectorIt.size(),1);
            for ( unsigned int i = 0; i < inputVectorIt.size(); i++ )
            {
                intensities( i,0 ) = inputVectorIt[ i ].Get();
            }

            float maxProbability = -1.0;
            int maxProbaIndex = -1;
            for ( unsigned int i = 0; i < m_GaussianModel.size(); i++ )
            {
                float proba = m_Alphas[i] * this->probability( intensities, m_GaussianModel[i] );
                if ( maxProbability < proba )
                {
                    maxProbability = proba;
                    maxProbaIndex = i+1;
                }
            }
            if ( maxProbaIndex > 0 )
            {
                classificationIt.Set(maxProbaIndex);
            }
            else
            {
                classificationIt.Set(m_GaussianModel.size()+1);
            }
        }
        ++classificationIt;
        ++maskIt;
        for ( unsigned int i = 0; i < m_ImagesVector.size(); i++ )
        {
            ++inputVectorIt[i];
        }
    }

}


template <typename TInput, typename TMask, typename TOutput>
double
ImageClassifierFilter<TInput, TMask, TOutput>::probability(DoubleVariableSizeMatrixType &point, GaussianFunctionType::Pointer model)
{
    DoubleVariableSizeMatrixType distance, distanceT, value, sigma, mu;
    sigma = model->GetCovariance();

    GaussianFunctionType::MeanVectorType muVect = model->GetMean();
    mu.SetSize(muVect.Size(),1);
    for(unsigned int j = 0; j < muVect.Size(); j++)
    {
        mu(j,0) = muVect[j];
    }

    double down = static_cast<double>(std::sqrt( std::pow(2*M_PI,muVect.Size()) * std::abs( vnl_determinant(sigma.GetVnlMatrix())) ));
    distance = point - mu;
    distanceT = distance.GetTranspose();
    value = distanceT * (sigma.GetInverse() * distance.GetVnlMatrix());

    if(value(0,0) != value(0,0))
        return 0.0; //NB index auparavant a 1, verif if ok

    double exponent =- 0.5* value(0,0);
    return exp(exponent)/down;
}

} //end of namespace anima
