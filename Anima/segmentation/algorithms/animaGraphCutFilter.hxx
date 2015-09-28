#pragma once

#include "animaGraphCutFilter.h"

namespace anima
{

template <typename TInput, typename TOutput>
void GraphCutFilter<TInput, TOutput>::SetInputImage1(const TInput* image)
{
    this->SetNthInput(0, const_cast<TInput*>(image));
    m_TLinksFilter->SetInputImage1( image );
    m_Graph3DFilter->SetInputImage1( image );
}

template <typename TInput, typename TOutput>
void GraphCutFilter<TInput, TOutput>::SetMask(const TMask* mask)
{
    this->SetNthInput(1, const_cast<TMask*>(mask));
    m_Graph3DFilter->SetMask( mask );
}

template <typename TInput, typename TOutput>
void GraphCutFilter<TInput, TOutput>::SetInputSeedSourcesMask(const TMask* mask)
{
    this->SetNthInput(m_NbInputs, const_cast<TMask*>(mask));
    m_IndexSourcesMask = m_NbInputs;
    m_NbInputs++;
    m_TLinksFilter->SetInputSeedSourcesMask( mask );
}

template <typename TInput, typename TOutput>
void GraphCutFilter<TInput, TOutput>::SetInputSeedSinksMask(const TMask* mask)
{
    this->SetNthInput(m_NbInputs, const_cast<TMask*>(mask));
    m_IndexSinksMask = m_NbInputs;
    m_NbInputs++;
    m_TLinksFilter->SetInputSeedSinksMask( mask );
}


template <typename TInput, typename TOutput>
void GraphCutFilter<TInput, TOutput>::SetInputSeedSourcesProba(const TSeedProba* mask)
{
    this->SetNthInput(m_NbInputs, const_cast<TSeedProba*>(mask));
    m_IndexSourcesProba = m_NbInputs;
    m_NbInputs++;
    m_TLinksFilter->SetInputSeedSourcesProba( mask );
}

template <typename TInput, typename TOutput>
void GraphCutFilter<TInput, TOutput>::SetInputSeedSinksProba(const TSeedProba* mask)
{
    this->SetNthInput(m_NbInputs, const_cast<TSeedProba*>(mask));
    m_IndexSinksProba = m_NbInputs;
    m_NbInputs++;
    m_TLinksFilter->SetInputSeedSinksProba( mask );
}

template <typename TInput, typename TOutput>
void GraphCutFilter<TInput, TOutput>::SetInputImage2(const TInput* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInput*>(image));
    m_IndexImage2 = m_NbInputs;
    m_NbInputs++;
    m_TLinksFilter->SetInputImage2( image );
    m_Graph3DFilter->SetInputImage2( image );
}

template <typename TInput, typename TOutput>
void GraphCutFilter<TInput, TOutput>::SetInputImage3(const TInput* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInput*>(image));
    m_IndexImage3 = m_NbInputs;
    m_NbInputs++;
    m_TLinksFilter->SetInputImage3( image );
    m_Graph3DFilter->SetInputImage3( image );
}
template <typename TInput, typename TOutput>
void GraphCutFilter<TInput, TOutput>::SetInputImage4(const TInput* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInput*>(image));
    m_IndexImage4 = m_NbInputs;
    m_NbInputs++;
    m_TLinksFilter->SetInputImage4( image );
    m_Graph3DFilter->SetInputImage4( image );
}

template <typename TInput, typename TOutput>
void GraphCutFilter<TInput, TOutput>::SetInputImage5(const TInput* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInput*>(image));
    m_IndexImage5 = m_NbInputs;
    m_NbInputs++;
    m_TLinksFilter->SetInputImage5( image );
    m_Graph3DFilter->SetInputImage5( image );
}



template <typename TInput, typename TOutput>
typename TInput::ConstPointer GraphCutFilter<TInput, TOutput>::GetInputImage1()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(0) );
}

template <typename TInput, typename TOutput>
itk::Image <unsigned char,3>::ConstPointer GraphCutFilter<TInput, TOutput>::GetMask()
{
    return static_cast< const TMask * >
            ( this->itk::ProcessObject::GetInput(1) );
}

template <typename TInput, typename TOutput>
itk::Image <unsigned char,3>::ConstPointer GraphCutFilter<TInput, TOutput>::GetInputSeedSourcesMask()
{
    return static_cast< const TMask * >
            ( this->itk::ProcessObject::GetInput(m_IndexSourcesMask) );
}

template <typename TInput, typename TOutput>
itk::Image <unsigned char,3>::ConstPointer GraphCutFilter<TInput, TOutput>::GetInputSeedSinksMask()
{
    return static_cast< const TMask * >
            ( this->itk::ProcessObject::GetInput(m_IndexSinksMask) );
}

template <typename TInput, typename TOutput>
itk::Image <double,3>::ConstPointer GraphCutFilter<TInput, TOutput>::GetInputSeedSourcesProba()
{
    return static_cast< const TSeedProba * >
            ( this->itk::ProcessObject::GetInput(m_IndexSourcesProba) );
}

template <typename TInput, typename TOutput>
itk::Image <double,3>::ConstPointer GraphCutFilter<TInput, TOutput>::GetInputSeedSinksProba()
{
    return static_cast< const TSeedProba * >
            ( this->itk::ProcessObject::GetInput(m_IndexSinksProba) );
}

template <typename TInput, typename TOutput>
typename TInput::ConstPointer GraphCutFilter<TInput, TOutput>::GetInputImage2()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage2) );
}

template <typename TInput, typename TOutput>
typename TInput::ConstPointer GraphCutFilter<TInput, TOutput>::GetInputImage3()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage3) );
}
template <typename TInput, typename TOutput>
typename TInput::ConstPointer GraphCutFilter<TInput, TOutput>::GetInputImage4()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage4) );
}

template <typename TInput, typename TOutput>
typename TInput::ConstPointer GraphCutFilter<TInput, TOutput>::GetInputImage5()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage5) );
}

template <typename TInputImage, typename TOutput>
itk::DataObject::Pointer GraphCutFilter <TInputImage, TOutput>::MakeOutput(unsigned int idx)
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

template <typename TInputImage, typename TOutput>
typename TOutput::Pointer GraphCutFilter <TInputImage, TOutput>::GetOutput()
{
    return dynamic_cast< TOutput* >( this->itk::ProcessObject::GetOutput(0) );
}
template <typename TInputImage, typename TOutput>
typename TOutput::Pointer GraphCutFilter <TInputImage, TOutput>::GetOutputBackground()
{
    return dynamic_cast< TOutput* >( this->itk::ProcessObject::GetOutput(1) );
}

template <typename TInput, typename TOutput>
void
GraphCutFilter<TInput, TOutput>
::WriteOutputs()
{
    if( m_OutputFilename != "" )
    {
        std::cout << "Writing graph cut output image to: " << m_OutputFilename << std::endl;
        anima::writeImage<TOutput>(m_OutputFilename, this->GetOutput());
    }
    if( m_OutputBackgroundFilename != "" )
    {
        std::cout << "Writing graph cut output image to: " << m_OutputBackgroundFilename << std::endl;
        anima::writeImage<TOutput>(m_OutputBackgroundFilename, this->GetOutputBackground());
    }
}

template <typename TInput, typename TOutput>
void
GraphCutFilter<TInput, TOutput>
::GenerateData()
{
    std::cout << "Computing graph cut..." << std::endl;

    m_TLinksFilter->SetAlpha( this->GetAlpha() );
    m_TLinksFilter->SetMultiVarSources( this->GetMultiVarSources() );
    m_TLinksFilter->SetMultiVarSinks( this->GetMultiVarSinks() );
    m_TLinksFilter->SetTLinkMode( this->GetTLinkMode() );
    m_TLinksFilter->SetVerbose( this->GetVerbose() );
    m_TLinksFilter->SetNumberOfThreads(this->GetNumberOfThreads());
    m_TLinksFilter->SetTol( m_Tol );

    m_Graph3DFilter->SetSigma( this->GetSigma() );
    m_Graph3DFilter->SetUseSpectralGradient( this->GetUseSpectralGradient() );
    m_Graph3DFilter->SetMatrix( this->GetMatrix() );
    m_Graph3DFilter->SetMatFilename( this->GetMatrixGradFilename() );
    m_Graph3DFilter->SetVerbose( this->GetVerbose() );
    m_Graph3DFilter->SetNumberOfThreads(this->GetNumberOfThreads());
    m_Graph3DFilter->SetTol( m_Tol );

    m_Graph3DFilter->GraftNthOutput( 0, this->GetOutput() );
    m_Graph3DFilter->GraftNthOutput( 1, this->GetOutputBackground() );

    try
    {
        m_Graph3DFilter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        exit(-1);
    }

    this->GraftNthOutput( 0 , m_Graph3DFilter->GetOutput() );
    this->GraftNthOutput( 1 , m_Graph3DFilter->GetOutputBackground() );
}


} //end of namespace anima
