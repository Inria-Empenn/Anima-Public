#pragma once

#include "animaComputeMahalanobisImagesFilter.h"

namespace anima
{

template <typename TInputImage, typename TMaskImage, typename TOutput>
void ComputeMahalanobisImagesFilter <TInputImage, TMaskImage, TOutput>::SetInputImage1(const TInputImage* image)
{
    this->SetNthInput(0, const_cast<TInputImage*>(image));
}

template <typename TInputImage, typename TMaskImage, typename TOutput>
void ComputeMahalanobisImagesFilter <TInputImage, TMaskImage, TOutput>::SetInputImage2(const TInputImage* image)
{
    this->SetNthInput(1, const_cast<TInputImage*>(image));
}

template <typename TInputImage, typename TMaskImage, typename TOutput>
void ComputeMahalanobisImagesFilter <TInputImage, TMaskImage, TOutput>::SetInputImage3(const TInputImage* image)
{
    this->SetNthInput(2, const_cast<TInputImage*>(image));
}

template <typename TInputImage, typename TMaskImage, typename TOutput>
void ComputeMahalanobisImagesFilter <TInputImage, TMaskImage, TOutput>::SetMask(const TMaskImage* mask)
{
    this->SetNthInput(3, const_cast<TMaskImage*>(mask));
}

template <typename TInputImage, typename TMaskImage, typename TOutput>
typename TInputImage::ConstPointer ComputeMahalanobisImagesFilter <TInputImage, TMaskImage, TOutput>::GetInputImage1()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(0) );
}

template <typename TInputImage, typename TMaskImage, typename TOutput>
typename TInputImage::ConstPointer ComputeMahalanobisImagesFilter <TInputImage, TMaskImage, TOutput>::GetInputImage2()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(1) );
}

template <typename TInputImage, typename TMaskImage, typename TOutput>
typename TInputImage::ConstPointer ComputeMahalanobisImagesFilter <TInputImage, TMaskImage, TOutput>::GetInputImage3()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(2) );
}

template <typename TInputImage, typename TMaskImage, typename TOutput>
typename TMaskImage::ConstPointer ComputeMahalanobisImagesFilter <TInputImage, TMaskImage, TOutput>::GetMask()
{
    return static_cast< const TMaskImage * >
            ( this->itk::ProcessObject::GetInput(3) );
}


template <typename TInputImage, typename TMaskImage, typename TOutput>
itk::DataObject::Pointer ComputeMahalanobisImagesFilter <TInputImage, TMaskImage, TOutput>::MakeOutput(unsigned int idx)
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
    case 2:
        output = ( TOutput::New() ).GetPointer();
        break;
    case 3:
        output = ( TOutput::New() ).GetPointer();
        break;
    case 4:
        output = ( TOutput::New() ).GetPointer();
        break;
    default:
        std::cerr << "No output " << idx << std::endl;
        output = NULL;
        break;
    }
    return output.GetPointer();
}

template <typename TInputImage, typename TMaskImage, typename TOutput>
typename TOutput::Pointer ComputeMahalanobisImagesFilter <TInputImage, TMaskImage, TOutput>::GetOutputMahaCSF()
{
    return dynamic_cast< OutputImageType * >( this->itk::ProcessObject::GetOutput(0) );
}
template <typename TInputImage, typename TMaskImage, typename TOutput>
typename TOutput::Pointer ComputeMahalanobisImagesFilter <TInputImage, TMaskImage, TOutput>::GetOutputMahaGM()
{
    return dynamic_cast< OutputImageType * >( this->itk::ProcessObject::GetOutput(1) );
}
template <typename TInputImage, typename TMaskImage, typename TOutput>
typename TOutput::Pointer ComputeMahalanobisImagesFilter <TInputImage, TMaskImage, TOutput>::GetOutputMahaWM()
{
    return dynamic_cast< OutputImageType * >( this->itk::ProcessObject::GetOutput(2) );
}
template <typename TInputImage, typename TMaskImage, typename TOutput>
typename TOutput::Pointer ComputeMahalanobisImagesFilter <TInputImage, TMaskImage, TOutput>::GetOutputMahaMinimum()
{
    return dynamic_cast< OutputImageType * >( this->itk::ProcessObject::GetOutput(3) );
}
template <typename TInputImage, typename TMaskImage, typename TOutput>
typename TOutput::Pointer ComputeMahalanobisImagesFilter <TInputImage, TMaskImage, TOutput>::GetOutputMahaMaximum()
{
    return dynamic_cast< OutputImageType * >( this->itk::ProcessObject::GetOutput(4) );
}

template <typename TInputImage, typename TMaskImage, typename TOutput>
void
ComputeMahalanobisImagesFilter <TInputImage, TMaskImage, TOutput>::WriteOutputs()
{
    if( m_OutputMahaCSFFilename != "" )
    {
        std::cout << "Writing mahalanobis CSF image to: " << m_OutputMahaCSFFilename << std::endl;
        anima::writeImage<TOutput>(m_OutputMahaCSFFilename, this->GetOutputMahaCSF());
    }
    if( m_OutputMahaGMFilename != "" )
    {
        std::cout << "Writing mahalanobis GM image to: " << m_OutputMahaGMFilename << std::endl;
        anima::writeImage<TOutput>(m_OutputMahaGMFilename, this->GetOutputMahaGM());
    }
    if( m_OutputMahaWMFilename != "" )
    {
        std::cout << "Writing mahalanobis WM image to: " << m_OutputMahaWMFilename << std::endl;
        anima::writeImage<TOutput>(m_OutputMahaWMFilename, this->GetOutputMahaWM());
    }
    if( m_OutputMahaMinimumFilename != "" )
    {
        std::cout << "Writing minimum mahalanobis image to: " << m_OutputMahaMinimumFilename << std::endl;
        anima::writeImage<TOutput>(m_OutputMahaMinimumFilename, this->GetOutputMahaMinimum());
    }
    if( m_OutputMahaMaximumFilename != "" )
    {
        std::cout << "Writing maximum mahalanobis image to: " << m_OutputMahaMaximumFilename << std::endl;
        anima::writeImage<TOutput>(m_OutputMahaMaximumFilename, this->GetOutputMahaMaximum());
    }
}

template <typename TInputImage, typename TMaskImage, typename TOutput>
void
ComputeMahalanobisImagesFilter <TInputImage, TMaskImage, TOutput>
::BeforeThreadedGenerateData()
{
    if(m_GaussianModel.size() != 3)
    {
        std::cout << "Error in mahalanobis filter: number of images do not correspond to solution size exiting..."<< std::endl;
        exit(-1);
    }

    m_InputImage_1_UC = ImageTypeUC::New();
    m_InputImage_2_UC = ImageTypeUC::New();
    m_InputImage_3_UC = ImageTypeUC::New();

    // Rescale images
    double desiredMinimum=0,desiredMaximum=255;

    typename RescaleFilterType::Pointer rescaleFilter1 = RescaleFilterType::New();
    rescaleFilter1->SetInput( this->GetInputImage1() );
    rescaleFilter1->SetOutputMinimum( desiredMinimum );
    rescaleFilter1->SetOutputMaximum( desiredMaximum );
    rescaleFilter1->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    rescaleFilter1->SetCoordinateTolerance(m_Tol);
    rescaleFilter1->SetDirectionTolerance(m_Tol);
    rescaleFilter1->UpdateLargestPossibleRegion();
    m_InputImage_1_UC = rescaleFilter1->GetOutput();

    typename RescaleFilterType::Pointer rescaleFilter2 = RescaleFilterType::New();
    rescaleFilter2->SetInput( this->GetInputImage2() );
    rescaleFilter2->SetOutputMinimum( desiredMinimum );
    rescaleFilter2->SetOutputMaximum( desiredMaximum );
    rescaleFilter2->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    rescaleFilter2->SetCoordinateTolerance(m_Tol);
    rescaleFilter2->SetDirectionTolerance(m_Tol);
    rescaleFilter2->UpdateLargestPossibleRegion();
    m_InputImage_2_UC = rescaleFilter2->GetOutput();

    typename RescaleFilterType::Pointer rescaleFilter3 = RescaleFilterType::New();
    rescaleFilter3->SetInput( this->GetInputImage3() );
    rescaleFilter3->SetOutputMinimum( desiredMinimum );
    rescaleFilter3->SetOutputMaximum( desiredMaximum );
    rescaleFilter3->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    rescaleFilter3->SetCoordinateTolerance(m_Tol);
    rescaleFilter3->SetDirectionTolerance(m_Tol);
    rescaleFilter3->UpdateLargestPossibleRegion();
    m_InputImage_3_UC = rescaleFilter3->GetOutput();
}


template <typename TInputImage, typename TMaskImage, typename TOutput>
void
ComputeMahalanobisImagesFilter <TInputImage, TMaskImage, TOutput>::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    const unsigned int                         Dimension = 3;
    typedef double                      PixelComponentType;
    typedef itk::RGBPixel<PixelComponentType>  PixelType;
    typedef itk::Image< PixelType, Dimension > ImageType;
    typedef itk::MahalanobisDistanceThresholdImageFunction< ImageType > FunctionType;
    typedef FunctionType::IndexType FunctionIndexType;

    OutputIteratorType mahaCSFIt(this->GetOutputMahaCSF(), outputRegionForThread );
    OutputIteratorType mahaGMIt(this->GetOutputMahaGM(), outputRegionForThread );
    OutputIteratorType mahaWMIt(this->GetOutputMahaWM(), outputRegionForThread );

    OutputIteratorType mahaMaxiIt(this->GetOutputMahaMaximum(), outputRegionForThread );
    OutputIteratorType mahaMiniIt(this->GetOutputMahaMinimum(), outputRegionForThread );

    ImageType::Pointer image = ImageType::New();
    image->SetRegions( this->GetMask()->GetLargestPossibleRegion() );
    image->CopyInformation( this->GetMask() );
    image->Allocate();
    itk::ImageRegionIterator<ImageType> itRGB (image, outputRegionForThread);

    ImageConstIteratorTypeUC It1(m_InputImage_1_UC, outputRegionForThread );
    ImageConstIteratorTypeUC It2(m_InputImage_2_UC, outputRegionForThread );
    ImageConstIteratorTypeUC It3(m_InputImage_3_UC, outputRegionForThread );

    while(!itRGB.IsAtEnd())
    {
        PixelType pix;
        pix[0] = It1.Get();
        pix[1] = It2.Get();
        pix[2] = It3.Get();
        itRGB.Set(pix);
        ++It1;
        ++It2;
        ++It3;
        ++itRGB;
    }

    FunctionType::Pointer function = FunctionType::New();
    FunctionType::CovarianceMatrixType Covariance( m_NbModalities, m_NbModalities );
    Covariance.fill( 0.0 );

    FunctionType::MeanVectorType Mean( m_NbModalities );
    function->SetInputImage( image );

    MaskConstIteratorType MaskImageIt ( this->GetMask(), outputRegionForThread );

    itk::Statistics::ChiSquareDistribution::Pointer chi = itk::Statistics::ChiSquareDistribution::New();
    itk::SizeValueType degree = (m_GaussianModel[0])->GetMean().Size();
    chi->SetDegreesOfFreedom(degree);

    while(!MaskImageIt.IsAtEnd())
    {
        mahaCSFIt.Set(0);
        mahaGMIt.Set(0);
        mahaWMIt.Set(0);
        mahaMaxiIt.Set(0);
        mahaMiniIt.Set(0);

        if(MaskImageIt.Get()!=0)
        {
            FunctionIndexType ind;
            ind[0] = MaskImageIt.GetIndex()[0];
            ind[1] = MaskImageIt.GetIndex()[1];
            ind[2] = MaskImageIt.GetIndex()[2];

            for(unsigned int i = 0; i < m_NbTissus; i++)
            {
                GaussianFunctionType::MeanVectorType mean = (m_GaussianModel[i])->GetMean();
                Mean[0] = mean[0];
                Mean[1] = mean[1];
                Mean[2] = mean[2];

                GaussianFunctionType::CovarianceMatrixType cov = (m_GaussianModel[i])->GetCovariance();
                Covariance[0][0] = cov[0][0]; Covariance[1][0] = cov[1][0]; Covariance[2][0] = cov[2][0];
                Covariance[0][1] = cov[0][1]; Covariance[1][1] = cov[1][1]; Covariance[2][1] = cov[2][1];
                Covariance[0][2] = cov[0][2]; Covariance[1][2] = cov[1][2]; Covariance[2][2] = cov[2][2];

                function->SetCovariance( Covariance );
                function->SetMean( Mean );
                double distanceD = function->EvaluateDistanceAtIndex( ind );
                double chiRes = chi->EvaluateCDF(distanceD);
                double distance = 1-chiRes;

                if(i==0)
                {
                    mahaCSFIt.Set(static_cast<OutputPixelType>(distance));
                }
                if(i==1)
                {
                    mahaGMIt.Set(static_cast<OutputPixelType>(distance));
                }
                if(i==2)
                {
                    mahaWMIt.Set(static_cast<OutputPixelType>(distance));
                }


                if(mahaMaxiIt.Get() < mahaCSFIt.Get())
                {
                    mahaMaxiIt.Set(mahaCSFIt.Get());
                }
                if(mahaMaxiIt.Get() < mahaGMIt.Get())
                {
                    mahaMaxiIt.Set(mahaGMIt.Get());
                }
                if(mahaMaxiIt.Get() < mahaWMIt.Get())
                {
                    mahaMaxiIt.Set(mahaWMIt.Get());
                }
                mahaMiniIt.Set(1.0 - mahaMaxiIt.Get());
            }
        }

        ++mahaCSFIt;
        ++mahaGMIt;
        ++mahaWMIt;
        ++mahaMaxiIt;
        ++mahaMiniIt;
        ++MaskImageIt;
    }
}

} //end of namespace anima
