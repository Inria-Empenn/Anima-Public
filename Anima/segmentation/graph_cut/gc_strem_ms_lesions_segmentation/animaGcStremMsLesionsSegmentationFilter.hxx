#include "animaGcStremMsLesionsSegmentationFilter.h"

namespace anima
{

template <typename TInputImage>
void GcStremMsLesionsSegmentationFilter <TInputImage>::SetMask(const TInputImage* mask)
{
    this->SetNthInput(m_NbInputs, const_cast<TInputImage*>(mask));
    m_IndexMask = m_NbInputs;
    m_NbInputs++;
}

template <typename TInputImage>
void GcStremMsLesionsSegmentationFilter <TInputImage>::SetInputImageT1(const TInputImage* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInputImage*>(image));
    m_IndexImageT1 = m_NbInputs;
    m_NbInputs++;
    m_Modalities++;
}

template <typename TInputImage>
void GcStremMsLesionsSegmentationFilter <TInputImage>::SetInputImageT2(const TInputImage* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInputImage*>(image));
    m_IndexImageT2 = m_NbInputs;
    m_NbInputs++;
    m_Modalities++;
}

template <typename TInputImage>
void GcStremMsLesionsSegmentationFilter <TInputImage>::SetInputImageDP(const TInputImage* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInputImage*>(image));
    m_IndexImageDP = m_NbInputs;
    m_NbInputs++;
    m_Modalities++;
}

template <typename TInputImage>
void GcStremMsLesionsSegmentationFilter <TInputImage>::SetInputImageFLAIR(const TInputImage* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInputImage*>(image));
    m_IndexImageFLAIR = m_NbInputs;
    m_NbInputs++;
    m_Modalities++;
}

template <typename TInputImage>
void GcStremMsLesionsSegmentationFilter <TInputImage>::SetInputImageT1Gd(const TInputImage* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInputImage*>(image));
    m_IndexImageT1Gd = m_NbInputs;
    m_NbInputs++;
    m_Modalities++;
}

template <typename TInputImage>
void GcStremMsLesionsSegmentationFilter <TInputImage>::SetInputCSFAtlas(const ImageTypeD* atlas)
{
    this->SetNthInput(m_NbInputs, const_cast<ImageTypeD*>(atlas));
    m_IndexAtlasCSF = m_NbInputs;
    m_NbInputs++;
}

template <typename TInputImage>
void GcStremMsLesionsSegmentationFilter <TInputImage>::SetInputGMAtlas(const ImageTypeD* atlas)
{
    this->SetNthInput(m_NbInputs, const_cast<ImageTypeD*>(atlas));
    m_IndexAtlasGM = m_NbInputs;
    m_NbInputs++;
}

template <typename TInputImage>
void GcStremMsLesionsSegmentationFilter <TInputImage>::SetInputWMAtlas(const ImageTypeD* atlas)
{
    this->SetNthInput(m_NbInputs, const_cast<ImageTypeD*>(atlas));
    m_IndexAtlasWM = m_NbInputs;
    m_NbInputs++;
}

template <typename TInputImage>
void GcStremMsLesionsSegmentationFilter <TInputImage>::SetInputLesionPrior(const ImageTypeD* image)
{
    this->SetNthInput(m_NbInputs, const_cast<ImageTypeD*>(image));
    m_IndexLesionPrior = m_NbInputs;
    m_NbInputs++;
}

template <typename TInputImage>
void GcStremMsLesionsSegmentationFilter <TInputImage>::SetSourcesMask(const ImageTypeUC* image)
{
    this->SetNthInput(m_NbInputs, const_cast<ImageTypeUC*>(image));
    m_IndexSourcesMask = m_NbInputs;
    m_NbInputs++;
}

template <typename TInputImage>
void GcStremMsLesionsSegmentationFilter <TInputImage>::SetSinksMask(const ImageTypeUC* image)
{
    this->SetNthInput(m_NbInputs, const_cast<ImageTypeUC*>(image));
    m_IndexSinksMask = m_NbInputs;
    m_NbInputs++;
}



template <typename TInputImage>
typename TInputImage::ConstPointer GcStremMsLesionsSegmentationFilter <TInputImage>::GetMask()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(m_IndexMask) );
}

template <typename TInputImage>
typename TInputImage::ConstPointer GcStremMsLesionsSegmentationFilter <TInputImage>::GetInputImageT1()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(m_IndexImageT1) );
}

template <typename TInputImage>
typename TInputImage::ConstPointer GcStremMsLesionsSegmentationFilter <TInputImage>::GetInputImageT2()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(m_IndexImageT2) );
}

template <typename TInputImage>
typename TInputImage::ConstPointer GcStremMsLesionsSegmentationFilter <TInputImage>::GetInputImageDP()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(m_IndexImageDP) );
}

template <typename TInputImage>
typename TInputImage::ConstPointer GcStremMsLesionsSegmentationFilter <TInputImage>::GetInputImageFLAIR()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(m_IndexImageFLAIR) );
}
template <typename TInputImage>
typename TInputImage::ConstPointer GcStremMsLesionsSegmentationFilter <TInputImage>::GetInputImageT1Gd()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(m_IndexImageT1Gd) );
}


template <typename TInputImage>
itk::Image <double,3>::ConstPointer GcStremMsLesionsSegmentationFilter <TInputImage>::GetInputCSFAtlas()
{
    return static_cast< const ImageTypeD * >
            ( this->itk::ProcessObject::GetInput(m_IndexAtlasCSF) );
}

template <typename TInputImage>
itk::Image <double,3>::ConstPointer GcStremMsLesionsSegmentationFilter <TInputImage>::GetInputGMAtlas()
{
    return static_cast< const ImageTypeD * >
            ( this->itk::ProcessObject::GetInput(m_IndexAtlasGM) );
}

template <typename TInputImage>
itk::Image <double,3>::ConstPointer GcStremMsLesionsSegmentationFilter <TInputImage>::GetInputWMAtlas()
{
    return static_cast< const ImageTypeD * >
            ( this->itk::ProcessObject::GetInput(m_IndexAtlasWM) );
}

template <typename TInputImage>
itk::Image <double,3>::ConstPointer GcStremMsLesionsSegmentationFilter <TInputImage>::GetInputLesionPrior()
{
    return static_cast< const ImageTypeD * > (this->itk::ProcessObject::GetInput(m_IndexLesionPrior));
}

template <typename TInputImage>
itk::Image <unsigned char,3>::ConstPointer GcStremMsLesionsSegmentationFilter <TInputImage>::GetSourcesMask()
{
    return static_cast< const ImageTypeUC * >
            ( this->itk::ProcessObject::GetInput(m_IndexSourcesMask) );
}
template <typename TInputImage>
itk::Image <unsigned char,3>::ConstPointer GcStremMsLesionsSegmentationFilter <TInputImage>::GetSinksMask()
{
    return static_cast< const ImageTypeUC * >
            ( this->itk::ProcessObject::GetInput(m_IndexSinksMask) );
}


template <typename TInputImage>
itk::DataObject::Pointer GcStremMsLesionsSegmentationFilter <TInputImage>::MakeOutput(unsigned int idx)
{
    itk::DataObject::Pointer output;

    if( (idx < 11) || (idx == 17) )
    {
        output = ( TOutputImage::New() ).GetPointer();
    }
    else if(( 11 <= idx ) && ( idx < 17 ))
    {
        output = ( ImageTypeD::New() ).GetPointer();
    }
    else
    {
        std::cerr << "No output " << idx << std::endl;
        output = NULL;
    }
    return output.GetPointer();
}

template <typename TInputImage>
itk::Image <unsigned char, 3>* GcStremMsLesionsSegmentationFilter <TInputImage>::GetOutputLesions()
{
    return dynamic_cast< TOutputImage* >( this->itk::ProcessObject::GetOutput(0) );
}
template <typename TInputImage>
itk::Image <unsigned char, 3>* GcStremMsLesionsSegmentationFilter <TInputImage>::GetOutputCSF()
{
    return dynamic_cast< TOutputImage* >( this->itk::ProcessObject::GetOutput(1) );
}
template <typename TInputImage>
itk::Image <unsigned char, 3>* GcStremMsLesionsSegmentationFilter <TInputImage>::GetOutputGM()
{
    return dynamic_cast< TOutputImage* >( this->itk::ProcessObject::GetOutput(2) );
}
template <typename TInputImage>
itk::Image <unsigned char, 3>* GcStremMsLesionsSegmentationFilter <TInputImage>::GetOutputWM()
{
    return dynamic_cast< TOutputImage* >( this->itk::ProcessObject::GetOutput(3) );
}
template <typename TInputImage>
itk::Image <unsigned char, 3>* GcStremMsLesionsSegmentationFilter <TInputImage>::GetOutputWholeSeg()
{
    return dynamic_cast< TOutputImage* >( this->itk::ProcessObject::GetOutput(4) );
}
template <typename TInputImage>
itk::Image <unsigned char, 3>* GcStremMsLesionsSegmentationFilter <TInputImage>::GetOutputStrem()
{
    return dynamic_cast< TOutputImage* >( this->itk::ProcessObject::GetOutput(5) );
}
template <typename TInputImage>
itk::Image <unsigned char, 3>* GcStremMsLesionsSegmentationFilter <TInputImage>::GetOutputStremCSF()
{
    return dynamic_cast< TOutputImage* >( this->itk::ProcessObject::GetOutput(6) );
}
template <typename TInputImage>
itk::Image <unsigned char, 3>* GcStremMsLesionsSegmentationFilter <TInputImage>::GetOutputStremGM()
{
    return dynamic_cast< TOutputImage* >( this->itk::ProcessObject::GetOutput(7) );
}
template <typename TInputImage>
itk::Image <unsigned char, 3>* GcStremMsLesionsSegmentationFilter <TInputImage>::GetOutputStremWM()
{
    return dynamic_cast< TOutputImage* >( this->itk::ProcessObject::GetOutput(8) );
}
template <typename TInputImage>
itk::Image <unsigned char, 3>* GcStremMsLesionsSegmentationFilter <TInputImage>::GetOutputIntensityImage1()
{
    return dynamic_cast< TOutputImage* >( this->itk::ProcessObject::GetOutput(9) );
}
template <typename TInputImage>
itk::Image <unsigned char, 3>* GcStremMsLesionsSegmentationFilter <TInputImage>::GetOutputIntensityImage2()
{
    return dynamic_cast< TOutputImage* >( this->itk::ProcessObject::GetOutput(10) );
}
template <typename TInputImage>
itk::Image <double,3>* GcStremMsLesionsSegmentationFilter <TInputImage>::GetOutputFuzzyObjectImage()
{
    return dynamic_cast< ImageTypeD* >( this->itk::ProcessObject::GetOutput(11) );
}
template <typename TInputImage>
itk::Image <double,3>* GcStremMsLesionsSegmentationFilter <TInputImage>::GetOutputMahaCSFImage()
{
    return dynamic_cast< ImageTypeD* >( this->itk::ProcessObject::GetOutput(12) );
}
template <typename TInputImage>
itk::Image <double,3>* GcStremMsLesionsSegmentationFilter <TInputImage>::GetOutputMahaGMImage()
{
    return dynamic_cast< ImageTypeD* >( this->itk::ProcessObject::GetOutput(13) );
}
template <typename TInputImage>
itk::Image <double,3>* GcStremMsLesionsSegmentationFilter <TInputImage>::GetOutputMahaWMImage()
{
    return dynamic_cast< ImageTypeD* >( this->itk::ProcessObject::GetOutput(14) );
}
template <typename TInputImage>
itk::Image <double,3>* GcStremMsLesionsSegmentationFilter <TInputImage>::GetOutputMahaMinimumImage()
{
    return dynamic_cast< ImageTypeD* >( this->itk::ProcessObject::GetOutput(15) );
}
template <typename TInputImage>
itk::Image <double,3>* GcStremMsLesionsSegmentationFilter <TInputImage>::GetOutputMahaMaximumImage()
{
    return dynamic_cast< ImageTypeD* >( this->itk::ProcessObject::GetOutput(16) );
}
template <typename TInputImage>
itk::Image <unsigned char, 3>* GcStremMsLesionsSegmentationFilter <TInputImage>::GetOutputGraphCut()
{
    return dynamic_cast< TOutputImage* >( this->itk::ProcessObject::GetOutput(17) );
}


template <typename TInputImage>
void
GcStremMsLesionsSegmentationFilter <TInputImage>::WriteOutputs()
{
    if( m_OutputLesionFilename != "" )
    {
        std::cout << "Writing Lesion output image to: " << m_OutputLesionFilename << std::endl;
        anima::writeImage<TOutputImage>(m_OutputLesionFilename, this->GetOutputLesions());
    }

    if( m_LesionSegmentationType!=manualGC )
    {
        m_MahalanobisFilter->WriteOutputs();
        if( m_OutputIntensityImage1Filename != "" )
        {
            std::cout << "Writing intensity output image 1 to: " << m_OutputIntensityImage1Filename << std::endl;
            anima::writeImage<TOutputImage>(m_OutputIntensityImage1Filename, this->GetOutputIntensityImage1());
        }
        if( m_OutputIntensityImage2Filename != "" )
        {
            std::cout << "Writing intensity output image 2 to: " << m_OutputIntensityImage2Filename << std::endl;
            anima::writeImage<TOutputImage>(m_OutputIntensityImage2Filename, this->GetOutputIntensityImage2());
        }
        if( m_OutputCSFFilename != "" )
        {
            std::cout << "Writing CSF output image to: " << m_OutputCSFFilename << std::endl;
            anima::writeImage<TOutputImage>(m_OutputCSFFilename, this->GetOutputCSF());
        }
        if( m_OutputGMFilename != "" )
        {
            std::cout << "Writing GM output image to: " << m_OutputGMFilename << std::endl;
            anima::writeImage<TOutputImage>(m_OutputGMFilename, this->GetOutputGM());
        }
        if( m_OutputWMFilename != "" )
        {
            std::cout << "Writing WM output image to: " << m_OutputWMFilename << std::endl;
            anima::writeImage<TOutputImage>(m_OutputWMFilename, this->GetOutputWM());
        }
        if( m_OutputWholeFilename != "" )
        {
            std::cout << "Writing segmentation output image to: " << m_OutputWholeFilename << std::endl;
            anima::writeImage<TOutputImage>(m_OutputWholeFilename, this->GetOutputWholeSeg());
        }
    }

    if( m_LesionSegmentationType==gcem || m_LesionSegmentationType==gcemAndManualGC )
    {
        if( m_OutputFuzzyObjectFilename != "" )
        {
            std::cout << "Writing fuzzy object output image to: " << m_OutputFuzzyObjectFilename << std::endl;
            anima::writeImage<ImageTypeD>(m_OutputFuzzyObjectFilename, this->GetOutputFuzzyObjectImage());
        }
    }

    if( m_LesionSegmentationType==strem )
    {
        if( m_OutputStremFilename != "" )
        {
            std::cout << "Writing strem output image to: " << m_OutputStremFilename << std::endl;
            anima::writeImage<TOutputImage>(m_OutputStremFilename, this->GetOutputStrem());
        }
        if( m_OutputStremCSFFilename != "" )
        {
            std::cout << "Writing CSF strem output image to: " << m_OutputStremCSFFilename << std::endl;
            anima::writeImage<TOutputImage>(m_OutputStremCSFFilename, this->GetOutputStremCSF());
        }
        if( m_OutputStremGMFilename != "" )
        {
            std::cout << "Writing GM strem output image to: " << m_OutputStremGMFilename << std::endl;
            anima::writeImage<TOutputImage>(m_OutputStremGMFilename, this->GetOutputStremGM());
        }
        if( m_OutputStremWMFilename != "" )
        {
            std::cout << "Writing WM strem output image to: " << m_OutputStremWMFilename << std::endl;
            anima::writeImage<TOutputImage>(m_OutputStremWMFilename, this->GetOutputStremWM());
        }
    }
    else
    {
        m_GraphCutFilter->WriteOutputs();
    }
}


template <typename TInputImage>
void
GcStremMsLesionsSegmentationFilter <TInputImage>::CheckInputImages()
{
    std::cout << "Check input entries..." << std::endl;
    if(m_LesionSegmentationType!=manualGC)
    {
        if(m_Modalities < 3)
            itkExceptionMacro("-- Error in anima::Gc_Strem_MS_Lesions_Segmentation: automatic segmentation requires 3 modalities, exiting...");

        if(this->GetInputImageT1().IsNull())
            itkExceptionMacro("-- Error in anima::Gc_Strem_MS_Lesions_Segmentation: automatic segmentation requires T1 image, exiting...");

        if(m_Modalities == 3)
        {
            m_UseT2 = false;
            m_UseDP = false;
            m_UseFLAIR = false;

            if(this->GetInputImageT2().IsNotNull()){m_UseT2=true;}
            if(this->GetInputImageDP().IsNotNull()){m_UseDP=true;}
            if(this->GetInputImageFLAIR().IsNotNull()){m_UseFLAIR=true;}
        }

        if(m_Modalities>3)
        {
            if( !( (m_UseT2 & m_UseDP & !m_UseFLAIR) || (m_UseT2 & !m_UseDP & m_UseFLAIR) || (!m_UseT2 & m_UseDP & m_UseFLAIR) ) )
                itkExceptionMacro("-- Error in anima::Gc_Strem_MS_Lesions_Segmentation: 2 images among T2, DP and FLAIR must be choosen for automatic segmentation, exiting...");
        }
    }
    else
    {
        if(this->GetSourcesMask().IsNull() || this->GetSinksMask().IsNull())
            itkExceptionMacro("-- Error in anima::Gc_Strem_MS_Lesions_Segmentation: seed mask missing, exiting...");
    }
}
template <typename TInputImage>
void
GcStremMsLesionsSegmentationFilter <TInputImage>::RescaleImages()
{
    std::cout << "Rescale input images..." << std::endl;

    ImageTypeUC::Pointer InputImage_T2_UC = ImageTypeUC::New();
    ImageTypeUC::Pointer InputImage_DP_UC = ImageTypeUC::New();
    ImageTypeUC::Pointer InputImage_FLAIR_UC = ImageTypeUC::New();

    double desiredMinimum=0,desiredMaximum=255;

    if(this->GetInputImageT1().IsNotNull())
    {
        typename RescaleFilterType::Pointer rescaleFilter1 = RescaleFilterType::New();
        rescaleFilter1->SetInput( this->GetInputImageT1() );
        rescaleFilter1->SetOutputMinimum( desiredMinimum );
        rescaleFilter1->SetOutputMaximum( desiredMaximum );
        rescaleFilter1->SetNumberOfThreads( this->GetNumberOfThreads() );
        rescaleFilter1->SetCoordinateTolerance( m_Tol );
        rescaleFilter1->SetDirectionTolerance( m_Tol );
        rescaleFilter1->UpdateLargestPossibleRegion();
        m_InputImage_T1_UC = rescaleFilter1->GetOutput();
    }
    if(this->GetInputImageT2().IsNotNull())
    {
        typename RescaleFilterType::Pointer rescaleFilter2 = RescaleFilterType::New();
        rescaleFilter2->SetInput( this->GetInputImageT2() );
        rescaleFilter2->SetOutputMinimum( desiredMinimum );
        rescaleFilter2->SetOutputMaximum( desiredMaximum );
        rescaleFilter2->SetNumberOfThreads( this->GetNumberOfThreads() );
        rescaleFilter2->SetCoordinateTolerance( m_Tol );
        rescaleFilter2->SetDirectionTolerance( m_Tol );
        rescaleFilter2->UpdateLargestPossibleRegion();
        InputImage_T2_UC = rescaleFilter2->GetOutput();
    }
    if(this->GetInputImageFLAIR().IsNotNull())
    {
        typename RescaleFilterType::Pointer rescaleFilter3 = RescaleFilterType::New();
        rescaleFilter3->SetInput( this->GetInputImageFLAIR() );
        rescaleFilter3->SetOutputMinimum( desiredMinimum );
        rescaleFilter3->SetOutputMaximum( desiredMaximum );
        rescaleFilter3->SetNumberOfThreads( this->GetNumberOfThreads() );
        rescaleFilter3->SetCoordinateTolerance( m_Tol );
        rescaleFilter3->SetDirectionTolerance( m_Tol );
        rescaleFilter3->UpdateLargestPossibleRegion();
        InputImage_FLAIR_UC = rescaleFilter3->GetOutput();
    }
    if(this->GetInputImageDP().IsNotNull())
    {
        typename RescaleFilterType::Pointer rescaleFilter4 = RescaleFilterType::New();
        rescaleFilter4->SetInput( this->GetInputImageDP() );
        rescaleFilter4->SetOutputMinimum( desiredMinimum );
        rescaleFilter4->SetOutputMaximum( desiredMaximum );
        rescaleFilter4->SetNumberOfThreads( this->GetNumberOfThreads() );
        rescaleFilter4->SetCoordinateTolerance( m_Tol );
        rescaleFilter4->SetDirectionTolerance( m_Tol );
        rescaleFilter4->UpdateLargestPossibleRegion();
        InputImage_DP_UC = rescaleFilter4->GetOutput();
    }

    // --------------- Create Brain Mask --------------- //
    m_MaskUC = ImageTypeUC::New();
    m_MaskUC->SetRegions( this->GetMask()->GetLargestPossibleRegion() );
    m_MaskUC->CopyInformation( this->GetMask() );
    m_MaskUC->Allocate();
    m_MaskUC->FillBuffer(0);

    ImageIteratorTypeUC itMaskUC( m_MaskUC, m_MaskUC->GetLargestPossibleRegion() );
    ImageIteratorTypeConstInput itInputImage_Mask( this->GetMask(), this->GetMask()->GetLargestPossibleRegion() );

    while (!itMaskUC.IsAtEnd())
    {
        if(itInputImage_Mask.Get()!=0)
        {
            itMaskUC.Set(1);
        }
        ++itMaskUC;
        ++itInputImage_Mask;
    }

    if(m_LesionSegmentationType!=manualGC)
    {
        if(m_UseT2 && m_UseDP)
        {
            m_InputImage_1_UC = InputImage_T2_UC;
            m_InputImage_2_UC = InputImage_DP_UC;
        }
        if(m_UseT2 && m_UseFLAIR)
        {
            m_InputImage_1_UC = InputImage_T2_UC;
            m_InputImage_2_UC = InputImage_FLAIR_UC;
        }
        if(m_UseDP && m_UseFLAIR)
        {
            m_InputImage_1_UC = InputImage_DP_UC;
            m_InputImage_2_UC = InputImage_FLAIR_UC;

            if(m_InitMethodType==1)
            {
                m_InputImage_1_UC = InputImage_FLAIR_UC;
                m_InputImage_2_UC = InputImage_DP_UC;
            }
        }
    }
}

template <typename TInputImage>
void
GcStremMsLesionsSegmentationFilter <TInputImage>::ComputeAutomaticInitialization()
{
    ComputeSolutionType::Pointer ComputeSolutionProcess = ComputeSolutionType::New();

    ComputeSolutionProcess->SetMask( m_MaskUC );
    ComputeSolutionProcess->SetInputImage1( m_InputImage_T1_UC );
    ComputeSolutionProcess->SetInputImage2( m_InputImage_1_UC );
    ComputeSolutionProcess->SetInputImage3( m_InputImage_2_UC );

    ComputeSolutionProcess->SetInputCSFAtlas( this->GetInputCSFAtlas() );
    ComputeSolutionProcess->SetInputGMAtlas( this->GetInputGMAtlas() );
    ComputeSolutionProcess->SetInputWMAtlas( this->GetInputWMAtlas() );

    ComputeSolutionProcess->SetSolutionWriteFilename( m_SolutionWriteFilename );
    ComputeSolutionProcess->SetSolutionReadFilename( m_SolutionReadFilename );
    ComputeSolutionProcess->SetGaussianModel( m_GaussianModel );
    ComputeSolutionProcess->SetAlphas( m_Alphas );

    ComputeSolutionProcess->SetUseT2( m_UseT2 );
    ComputeSolutionProcess->SetUseDP( m_UseDP );
    ComputeSolutionProcess->SetUseFLAIR( m_UseFLAIR );

    ComputeSolutionProcess->SetInitMethodType( m_InitMethodType );
    ComputeSolutionProcess->SetRejRatioHierar( m_RejRatioHierar );
    ComputeSolutionProcess->SetMinDistance( m_MinDistance );
    ComputeSolutionProcess->SetEmIter( m_EmIter );
    ComputeSolutionProcess->SetRejRatio( m_RejRatio );
    ComputeSolutionProcess->SetEmIter_concentration( m_EmIter_concentration );
    ComputeSolutionProcess->SetEM_before_concentration( m_EM_before_concentration );

    ComputeSolutionProcess->SetVerbose( m_Verbose );

    try
    {
        ComputeSolutionProcess->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        exit(-1);
    }

    m_GaussianModel = ComputeSolutionProcess->GetGaussianModel();

    std::cout << "Computing mahalanobis images..." << std::endl;

    m_MahalanobisFilter->SetGaussianModel( ComputeSolutionProcess->GetGaussianModel() );
    m_MahalanobisFilter->SetMask( m_MaskUC );
    m_MahalanobisFilter->SetInputImage1( m_InputImage_T1_UC );
    m_MahalanobisFilter->SetInputImage2( m_InputImage_1_UC );
    m_MahalanobisFilter->SetInputImage3( m_InputImage_2_UC );
    m_MahalanobisFilter->SetTol( m_Tol );
    m_MahalanobisFilter->SetNumberOfThreads( this->GetNumberOfThreads() );

    m_MahalanobisFilter->SetOutputMahaCSFFilename( m_OutputMahaCSFFilename );
    m_MahalanobisFilter->SetOutputMahaGMFilename( m_OutputMahaGMFilename );
    m_MahalanobisFilter->SetOutputMahaWMFilename( m_OutputMahaWMFilename );
    m_MahalanobisFilter->SetOutputMahaMaximumFilename( m_OutputMahaMaximumFilename );
    m_MahalanobisFilter->SetOutputMahaMinimumFilename ( m_OutputMahaMinimumFilename );

    m_MahalanobisFilter->GraftNthOutput( 0, this->GetOutputMahaCSFImage() );
    m_MahalanobisFilter->GraftNthOutput( 1, this->GetOutputMahaGMImage() );
    m_MahalanobisFilter->GraftNthOutput( 2, this->GetOutputMahaWMImage() );
    m_MahalanobisFilter->GraftNthOutput( 3, this->GetOutputMahaMinimumImage() );
    m_MahalanobisFilter->GraftNthOutput( 4, this->GetOutputMahaMaximumImage() );
    m_MahalanobisFilter->Update();
    this->GraftNthOutput( 12 , m_MahalanobisFilter->GetOutputMahaCSF() );
    this->GraftNthOutput( 13 , m_MahalanobisFilter->GetOutputMahaGM() );
    this->GraftNthOutput( 14 , m_MahalanobisFilter->GetOutputMahaWM() );
    this->GraftNthOutput( 15 , m_MahalanobisFilter->GetOutputMahaMinimum() );
    this->GraftNthOutput( 16 , m_MahalanobisFilter->GetOutputMahaMaximum() );

}

template <typename TInputImage>
void
GcStremMsLesionsSegmentationFilter <TInputImage>::StremThreshold()
{
    std::cout << "Compute strem threshold images..." << std::endl;

    const unsigned char insideValue = 1;
    const unsigned char outsideValue = 0;

    BinaryThresholdImageFilterType_F_UC::Pointer thresholdFilterCSF  = BinaryThresholdImageFilterType_F_UC::New();
    thresholdFilterCSF->SetInput( m_MahalanobisFilter->GetOutputMahaCSF() );
    thresholdFilterCSF->SetUpperThreshold( m_MahalanobisThCSF );
    thresholdFilterCSF->SetInsideValue( insideValue );
    thresholdFilterCSF->SetOutsideValue( outsideValue );
    thresholdFilterCSF->SetNumberOfThreads( this->GetNumberOfThreads() );
    thresholdFilterCSF->SetCoordinateTolerance( m_Tol );
    thresholdFilterCSF->SetDirectionTolerance( m_Tol );

    MaskFilterType_UC_UC::Pointer maskFilterCSF = MaskFilterType_UC_UC::New();
    maskFilterCSF->SetInput( thresholdFilterCSF->GetOutput()) ;
    maskFilterCSF->SetMaskImage( m_MaskUC );
    maskFilterCSF->SetNumberOfThreads( this->GetNumberOfThreads() );
    maskFilterCSF->SetCoordinateTolerance( m_Tol );
    maskFilterCSF->SetDirectionTolerance( m_Tol );
    maskFilterCSF->GraftOutput( this->GetOutputStremCSF() );
    maskFilterCSF->Update();
    this->GraftNthOutput( 6 , maskFilterCSF->GetOutput() );

    BinaryThresholdImageFilterType_F_UC::Pointer thresholdFilterGM  = BinaryThresholdImageFilterType_F_UC::New();
    thresholdFilterGM->SetInput( m_MahalanobisFilter->GetOutputMahaGM() );
    thresholdFilterGM->SetUpperThreshold( m_MahalanobisThGM );
    thresholdFilterGM->SetInsideValue( insideValue );
    thresholdFilterGM->SetOutsideValue( outsideValue );
    thresholdFilterGM->SetNumberOfThreads( this->GetNumberOfThreads() );
    thresholdFilterGM->SetCoordinateTolerance( m_Tol );
    thresholdFilterGM->SetDirectionTolerance( m_Tol );

    MaskFilterType_UC_UC::Pointer maskFilterGM = MaskFilterType_UC_UC::New();
    maskFilterGM->SetInput( thresholdFilterGM->GetOutput() );
    maskFilterGM->SetMaskImage( m_MaskUC );
    maskFilterGM->SetNumberOfThreads( this->GetNumberOfThreads() );
    maskFilterGM->SetCoordinateTolerance( m_Tol );
    maskFilterGM->SetDirectionTolerance( m_Tol );
    maskFilterGM->GraftOutput( this->GetOutputStremGM() );
    maskFilterGM->Update();
    this->GraftNthOutput( 7 , maskFilterGM->GetOutput() );

    BinaryThresholdImageFilterType_F_UC::Pointer thresholdFilterWM  = BinaryThresholdImageFilterType_F_UC::New();
    thresholdFilterWM->SetInput( m_MahalanobisFilter->GetOutputMahaWM() );
    thresholdFilterWM->SetUpperThreshold( m_MahalanobisThWM );
    thresholdFilterWM->SetInsideValue( insideValue );
    thresholdFilterWM->SetOutsideValue( outsideValue );
    thresholdFilterWM->SetNumberOfThreads( this->GetNumberOfThreads() );
    thresholdFilterWM->SetCoordinateTolerance( m_Tol );
    thresholdFilterWM->SetDirectionTolerance( m_Tol );

    MaskFilterType_UC_UC::Pointer maskFilterWM = MaskFilterType_UC_UC::New();
    maskFilterWM->SetInput( thresholdFilterWM->GetOutput() );
    maskFilterWM->SetMaskImage( m_MaskUC );
    maskFilterWM->SetNumberOfThreads( this->GetNumberOfThreads() );
    maskFilterWM->SetCoordinateTolerance( m_Tol );
    maskFilterWM->SetDirectionTolerance( m_Tol );
    maskFilterWM->GraftOutput( this->GetOutputStremWM() );
    maskFilterWM->Update();
    this->GraftNthOutput( 8 , maskFilterWM->GetOutput() );

    MinimumFilterTypeUC::Pointer filtermin1 = MinimumFilterTypeUC::New();
    MinimumFilterTypeUC::Pointer filtermin2 = MinimumFilterTypeUC::New();

    filtermin1->SetInput1( maskFilterCSF->GetOutput() );
    filtermin1->SetInput2( maskFilterGM->GetOutput() );
    filtermin1->SetNumberOfThreads( this->GetNumberOfThreads() );
    filtermin1->SetCoordinateTolerance( m_Tol);
    filtermin1->SetDirectionTolerance( m_Tol);

    filtermin2->SetInput1( filtermin1->GetOutput() );
    filtermin2->SetInput2( maskFilterWM->GetOutput() );
    filtermin2->SetNumberOfThreads( this->GetNumberOfThreads() );
    filtermin2->SetCoordinateTolerance( m_Tol );
    filtermin2->SetDirectionTolerance( m_Tol );
    filtermin2->GraftOutput( this->GetOutputStrem() );
    filtermin2->Update();
    this->GraftNthOutput( 5 , filtermin2->GetOutput() );
}

template <typename TInputImage>
void
GcStremMsLesionsSegmentationFilter <TInputImage>::CreateFuzzyRuleImage()
{
    std::cout << "Computing fuzzy rules..." << std::endl;

    MinimumFilterTypeF::Pointer filtermin1 = MinimumFilterTypeF::New();
    MinimumFilterTypeF::Pointer filtermin2 = MinimumFilterTypeF::New();

    const double minimumOutputValue = 0;
    const double maximumOutputValue = 1;

    const unsigned int indexWM = 2;
    GaussianFunctionType::MeanVectorType mean = (m_GaussianModel[indexWM])->GetMean();
    GaussianFunctionType::CovarianceMatrixType covar = (m_GaussianModel[indexWM])->GetCovariance();

    const unsigned int indexImage1 = 1;
    unsigned char intensityWindowLowerInputValue1 = static_cast<unsigned char>( mean[indexImage1] + m_FuzzyRuleMin * std::sqrt(covar[indexImage1][indexImage1]) );
    unsigned char intensityWindowUpperInputValue1 = static_cast<unsigned char>( mean[indexImage1] + m_FuzzyRuleMax * std::sqrt(covar[indexImage1][indexImage1]) );

    IntensityWindowingImageFilterType::Pointer intensityWindowFilter1 = IntensityWindowingImageFilterType::New();
    intensityWindowFilter1->SetInput( m_InputImage_1_UC );
    intensityWindowFilter1->SetWindowMinimum( intensityWindowLowerInputValue1 );
    intensityWindowFilter1->SetWindowMaximum( intensityWindowUpperInputValue1 );
    intensityWindowFilter1->SetOutputMinimum( minimumOutputValue );
    intensityWindowFilter1->SetOutputMaximum( maximumOutputValue );
    intensityWindowFilter1->SetNumberOfThreads( this->GetNumberOfThreads() );
    intensityWindowFilter1->SetCoordinateTolerance( m_Tol );
    intensityWindowFilter1->SetDirectionTolerance( m_Tol );

    const unsigned int indexImage2 = 2;
    unsigned char intensityWindowLowerInputValue2 = static_cast<unsigned char>( mean[indexImage2] + m_FuzzyRuleMin * std::sqrt(covar[indexImage2][indexImage2]) );
    unsigned char intensityWindowUpperInputValue2 = static_cast<unsigned char>( mean[indexImage2] + m_FuzzyRuleMax * std::sqrt(covar[indexImage2][indexImage2]) );

    IntensityWindowingImageFilterType::Pointer intensityWindowFilter2 = IntensityWindowingImageFilterType::New();
    intensityWindowFilter2->SetInput( m_InputImage_2_UC );
    intensityWindowFilter2->SetWindowMinimum( intensityWindowLowerInputValue2 );
    intensityWindowFilter2->SetWindowMaximum( intensityWindowUpperInputValue2 );
    intensityWindowFilter2->SetOutputMinimum( minimumOutputValue );
    intensityWindowFilter2->SetOutputMaximum( maximumOutputValue );
    intensityWindowFilter2->SetNumberOfThreads( this->GetNumberOfThreads() );
    intensityWindowFilter2->SetCoordinateTolerance( m_Tol );
    intensityWindowFilter2->SetDirectionTolerance( m_Tol );

    filtermin1->SetInput1( intensityWindowFilter1->GetOutput() );
    filtermin1->SetInput2( intensityWindowFilter2->GetOutput() );
    filtermin1->SetNumberOfThreads( this->GetNumberOfThreads() );
    filtermin1->SetCoordinateTolerance( m_Tol );
    filtermin1->SetDirectionTolerance( m_Tol );

    typename ImageTypeD::Pointer mahaMinimum = ImageTypeD::New();
    mahaMinimum->Graft(m_MahalanobisFilter->GetOutputMahaMinimum());

    if (this->GetInputLesionPrior().IsNotNull())
    {
        typedef itk::ImageRegionIterator <ImageTypeD> IteratorType;
        IteratorType mahaMinimumIterator(m_MahalanobisFilter->GetOutputMahaMinimum(),m_MahalanobisFilter->GetOutputMahaMinimum()->GetLargestPossibleRegion());
        typedef itk::ImageRegionConstIterator <ImageTypeD> ConstIteratorType;
        ConstIteratorType lesionPriorIterator(this->GetInputLesionPrior(),m_MahalanobisFilter->GetOutputMahaMinimum()->GetLargestPossibleRegion());

        while (!lesionPriorIterator.IsAtEnd())
        {
            double mahaValue = mahaMinimumIterator.Get();
            double lesionPriorValue = lesionPriorIterator.Get();

            mahaMinimumIterator.Set(std::sqrt(mahaValue * lesionPriorValue));

            ++mahaMinimumIterator;
            ++lesionPriorIterator;
        }
    }

    filtermin2->SetInput1( filtermin1->GetOutput() );
    filtermin2->SetInput2( mahaMinimum );
    filtermin2->SetNumberOfThreads( this->GetNumberOfThreads() );
    filtermin2->SetCoordinateTolerance( m_Tol );
    filtermin2->SetDirectionTolerance( m_Tol );

    MaskFilterType_F_UC::Pointer maskFilter = MaskFilterType_F_UC::New();
    maskFilter->SetInput( filtermin2->GetOutput() ) ;
    maskFilter->SetMaskImage( m_MaskUC );
    maskFilter->SetNumberOfThreads( this->GetNumberOfThreads() );
    maskFilter->SetCoordinateTolerance( m_Tol );
    maskFilter->SetDirectionTolerance( m_Tol );
    maskFilter->GraftOutput( this->GetOutputFuzzyObjectImage() );
    maskFilter->Update();
    this->GraftNthOutput( 11 , maskFilter->GetOutput() );

    m_FuzzyObject = ImageTypeD::New();
    m_FuzzyObject->Graft(maskFilter->GetOutput());
}

template <typename TInputImage>
void
GcStremMsLesionsSegmentationFilter <TInputImage>::GraphCut()
{
    m_GraphCutFilter->SetUseSpectralGradient( m_UseSpecGrad );
    m_GraphCutFilter->SetAlpha( m_Alpha );
    m_GraphCutFilter->SetSigma( m_Sigma );
    m_GraphCutFilter->SetMultiVarSources( m_MultiVarSources );
    m_GraphCutFilter->SetMultiVarSinks( m_MultiVarSinks );
    m_GraphCutFilter->SetMatrixGradFilename( m_MatrixGradFilename );
    m_GraphCutFilter->SetOutputFilename( m_OutputGCFilename );
    m_GraphCutFilter->SetVerbose( m_Verbose );
    m_GraphCutFilter->SetNumberOfThreads(this->GetNumberOfThreads());
    m_GraphCutFilter->SetTol( m_Tol );

    unsigned int index = 0;
    if(this->GetInputImageT1().IsNotNull()){ m_GraphCutFilter->SetInputImage(index, this->GetInputImageT1() ); ++index;}
    if(this->GetInputImageT2().IsNotNull()){ m_GraphCutFilter->SetInputImage(index, this->GetInputImageT2() ); ++index;}
    if(this->GetInputImageDP().IsNotNull()){ m_GraphCutFilter->SetInputImage(index, this->GetInputImageDP() ); ++index;}
    if(this->GetInputImageFLAIR().IsNotNull()){ m_GraphCutFilter->SetInputImage(index, this->GetInputImageFLAIR() ); ++index;}
    if(this->GetInputImageT1Gd().IsNotNull()){ m_GraphCutFilter->SetInputImage(index, this->GetInputImageT1Gd() ); ++index;}

    m_GraphCutFilter->SetMask( m_MaskUC );

    if( m_LesionSegmentationType==manualGC )
    {
        m_GraphCutFilter->SetTLinkMode( singleGaussianTLink );
        m_GraphCutFilter->SetInputSeedSourcesMask( this->GetSourcesMask() );
        m_GraphCutFilter->SetInputSeedSinksMask( this->GetSinksMask() );
    }
    else
    {
        this->CreateFuzzyRuleImage();

        // Compute TLinks if necessary
        m_TLinksFilter->SetAlpha( m_Alpha );
        m_TLinksFilter->SetMultiVarSources( m_MultiVarSources );
        m_TLinksFilter->SetMultiVarSinks( m_MultiVarSinks );
        m_TLinksFilter->SetTLinkMode( singleGaussianTLink );
        m_TLinksFilter->SetNumberOfThreads( this->GetNumberOfThreads() );
        m_TLinksFilter->SetTol( m_Tol );

        index = 0;
        if(this->GetInputImageT1().IsNotNull()){m_TLinksFilter->SetInputImage(index, this->GetInputImageT1() ); ++index;}
        if(this->GetInputImageT2().IsNotNull()){m_TLinksFilter->SetInputImage(index, this->GetInputImageT2() ); ++index;}
        if(this->GetInputImageDP().IsNotNull()){m_TLinksFilter->SetInputImage(index, this->GetInputImageDP() ); ++index;}
        if(this->GetInputImageFLAIR().IsNotNull()){m_TLinksFilter->SetInputImage(index, this->GetInputImageFLAIR() ); ++index;}
        if(this->GetInputImageT1Gd().IsNotNull()){m_TLinksFilter->SetInputImage(index, this->GetInputImageT1Gd() ); ++index;}

        m_TLinksFilter->SetInputSeedSourcesMask( this->GetSourcesMask() );
        m_TLinksFilter->SetInputSeedSinksMask( this->GetSinksMask() );

        if(this->GetSourcesMask().IsNotNull() && m_LesionSegmentationType==gcemAndManualGC)
        {
            m_FilterMaxSources->SetInput1( m_FuzzyObject );
            m_FilterMaxSources->SetInput2( m_TLinksFilter->GetOutputSources() );
            m_FilterMaxSources->SetNumberOfThreads( this->GetNumberOfThreads() );
            m_FilterMaxSources->SetCoordinateTolerance( m_Tol );
            m_FilterMaxSources->SetDirectionTolerance( m_Tol );
        }
        if(this->GetSinksMask().IsNotNull() && m_LesionSegmentationType==gcemAndManualGC)
        {
            m_FilterMaxSinks->SetInput1( m_MahalanobisFilter->GetOutputMahaMaximum() );
            m_FilterMaxSinks->SetInput2( m_TLinksFilter->GetOutputSinks() );
            m_FilterMaxSinks->SetNumberOfThreads( this->GetNumberOfThreads() );
            m_FilterMaxSinks->SetCoordinateTolerance( m_Tol );
            m_FilterMaxSinks->SetDirectionTolerance( m_Tol );
        }

        m_GraphCutFilter->SetTLinkMode( stremTLink );
        if(this->GetSourcesMask().IsNotNull() && m_LesionSegmentationType==gcemAndManualGC)
        {
            m_GraphCutFilter->SetInputSeedSourcesProba( m_FilterMaxSources->GetOutput() );
        }
        else
        {
            m_GraphCutFilter->SetInputSeedSourcesProba( m_FuzzyObject );
        }

        if(this->GetSinksMask().IsNotNull() && m_LesionSegmentationType==gcemAndManualGC)
        {
            m_GraphCutFilter->SetInputSeedSinksProba( m_FilterMaxSinks->GetOutput() );
        }
        else if (this->GetInputLesionPrior().IsNotNull())
        {
            typename ImageTypeD::Pointer sinkProbas = ImageTypeD::New();
            sinkProbas->Graft(m_MahalanobisFilter->GetOutputMahaMaximum());

            typedef itk::ImageRegionIterator <ImageTypeD> IteratorType;
            IteratorType sinkProbasIterator(m_MahalanobisFilter->GetOutputMahaMaximum(),m_MahalanobisFilter->GetOutputMahaMaximum()->GetLargestPossibleRegion());
            typedef itk::ImageRegionConstIterator <ImageTypeD> ConstIteratorType;
            ConstIteratorType lesionPriorIterator(this->GetInputLesionPrior(),m_MahalanobisFilter->GetOutputMahaMaximum()->GetLargestPossibleRegion());

            while (!lesionPriorIterator.IsAtEnd())
            {
                double mahaValue = sinkProbasIterator.Get();
                double lesionPriorValue = lesionPriorIterator.Get();

                sinkProbasIterator.Set(std::sqrt(mahaValue * (1.0 - lesionPriorValue)));

                ++sinkProbasIterator;
                ++lesionPriorIterator;
            }

            m_GraphCutFilter->SetInputSeedSinksProba(sinkProbas);
        }
        else
        {
            m_GraphCutFilter->SetInputSeedSinksProba(m_MahalanobisFilter->GetOutputMahaMaximum());
        }
    }
    m_GraphCutFilter->GraftOutput( this->GetOutputGraphCut() );
    m_GraphCutFilter->Update();
    this->GraftNthOutput( 17 , m_GraphCutFilter->GetOutput() );
}

template <typename TInputImage>
void
GcStremMsLesionsSegmentationFilter <TInputImage>::ApplyHeuristicRules()
{
    std::cout << "Computing heuristic rules..." << std::endl;

    ImageTypeUC::Pointer OutliersIntense = ImageTypeUC::New();
    unsigned int IndexT2ImageInModel, IndexDPImageInModel, IndexFLAIRImageInModel;
    const unsigned char insideValue = 1;
    const unsigned char outsideValue = 0;

    if( m_LesionSegmentationType != manualGC )
    {
        std::cout << "Computing hyper-intensity images..." << std::endl;

        /**
         * Intensity Rule: MS lesions are known to be hyperintense compared to the WM intensity on T2-w and PD-w and FLAIR sequences.
         * We use the information given by the NABT model to define hyper-intensity.
         * A voxel is considered to be hyper-intense for a given sequence if its intensity y is greater than a threshold T
         * that is defined by the probability of the Gaussian distribution and an intensity factor intensityFactor:
         * For example the hyperintensity threshold T in T2 is defined as:
         * T = meanT2 + intensityFactor * stdT2
         * meanT2 and stdT2 being the mean and standard deviation of the white matter in the sequence T2.
         * If the voxel is not considered hyper-intense on T2-w, PD and FLAIR, it is discarded as a lesion voxel.
         */

        GaussianFunctionType::CovarianceMatrixType covarWM = (m_GaussianModel[m_IndexWMinModel])->GetCovariance();
        GaussianFunctionType::MeanVectorType meanWM = (m_GaussianModel[m_IndexWMinModel])->GetMean();

        if(m_UseT2 && m_UseDP)
        {
            IndexT2ImageInModel = 1;
            double wmMeanT2 = meanWM [IndexT2ImageInModel];
            double wmCovarT2 = covarWM( IndexT2ImageInModel, IndexT2ImageInModel );
            m_HyperIntensityThreshold1 = wmMeanT2 + m_IntensityT2Factor * sqrt( wmCovarT2 );

            IndexDPImageInModel = 2;
            double wmMeanDP = meanWM[IndexDPImageInModel];
            double wmCovarDP = covarWM( IndexDPImageInModel, IndexDPImageInModel );
            m_HyperIntensityThreshold2 = wmMeanDP + m_IntensityDPFactor * sqrt( wmCovarDP );
        }
        if(m_UseT2 && m_UseFLAIR)
        {
            IndexT2ImageInModel = 1;
            double wmMeanT2 = meanWM[IndexT2ImageInModel];
            double wmCovarT2 = covarWM( IndexT2ImageInModel, IndexT2ImageInModel );
            m_HyperIntensityThreshold1 = wmMeanT2 + m_IntensityT2Factor * sqrt( wmCovarT2 );

            IndexFLAIRImageInModel = 2;
            double wmMeanFLAIR = meanWM[IndexFLAIRImageInModel];
            double wmCovarFLAIR = covarWM( IndexFLAIRImageInModel, IndexFLAIRImageInModel );
            m_HyperIntensityThreshold2 = wmMeanFLAIR + m_IntensityFLAIRFactor * sqrt( wmCovarFLAIR );
        }
        if(m_UseDP && m_UseFLAIR)
        {
            IndexDPImageInModel = 1;
            double wmMeanDP = meanWM[IndexDPImageInModel];
            double wmCovarDP = covarWM( IndexDPImageInModel, IndexDPImageInModel );
            m_HyperIntensityThreshold1 = wmMeanDP + m_IntensityDPFactor * sqrt( wmCovarDP );

            IndexFLAIRImageInModel = 2;
            double wmMeanFLAIR = meanWM[IndexFLAIRImageInModel];
            double wmCovarFLAIR = covarWM( IndexFLAIRImageInModel, IndexFLAIRImageInModel );
            m_HyperIntensityThreshold2 = wmMeanFLAIR + m_IntensityFLAIRFactor * sqrt( wmCovarFLAIR );

            if( m_InitMethodType == HierarchicalDP )
            {
                m_HyperIntensityThreshold1 = wmMeanFLAIR + m_IntensityFLAIRFactor * sqrt( wmCovarFLAIR );
                m_HyperIntensityThreshold2 = wmMeanDP + m_IntensityFLAIRFactor * sqrt( wmCovarDP );
            }
        }

        BinaryThresholdImageFilterType_UC_UC::Pointer IntensityFilter1  = BinaryThresholdImageFilterType_UC_UC::New();
        IntensityFilter1->SetInput( m_InputImage_1_UC );
        IntensityFilter1->SetLowerThreshold( m_HyperIntensityThreshold1 );
        IntensityFilter1->SetInsideValue( insideValue );
        IntensityFilter1->SetOutsideValue( outsideValue );
        IntensityFilter1->SetNumberOfThreads( this->GetNumberOfThreads() );
        IntensityFilter1->SetCoordinateTolerance( m_Tol );
        IntensityFilter1->SetDirectionTolerance( m_Tol );

        MaskFilterType_UC_UC::Pointer maskFilterIntensity1 = MaskFilterType_UC_UC::New();
        maskFilterIntensity1->SetInput( IntensityFilter1->GetOutput() ) ;
        maskFilterIntensity1->SetMaskImage( m_MaskUC );
        maskFilterIntensity1->SetNumberOfThreads( this->GetNumberOfThreads() );
        maskFilterIntensity1->SetCoordinateTolerance( m_Tol );
        maskFilterIntensity1->SetDirectionTolerance( m_Tol );
        maskFilterIntensity1->GraftOutput( this->GetOutputIntensityImage1() );
        maskFilterIntensity1->Update();
        this->GraftNthOutput( 9 , maskFilterIntensity1->GetOutput() );

        BinaryThresholdImageFilterType_UC_UC::Pointer IntensityFilter2 = BinaryThresholdImageFilterType_UC_UC::New();
        IntensityFilter2->SetInput( m_InputImage_2_UC );
        IntensityFilter2->SetLowerThreshold( m_HyperIntensityThreshold2 );
        IntensityFilter2->SetInsideValue( insideValue );
        IntensityFilter2->SetOutsideValue( outsideValue );
        IntensityFilter2->SetNumberOfThreads( this->GetNumberOfThreads() );
        IntensityFilter2->SetCoordinateTolerance( m_Tol );
        IntensityFilter2->SetDirectionTolerance( m_Tol );

        MaskFilterType_UC_UC::Pointer maskFilterIntensity2 = MaskFilterType_UC_UC::New();
        maskFilterIntensity2->SetInput( IntensityFilter2->GetOutput() ) ;
        maskFilterIntensity2->SetMaskImage( m_MaskUC );
        maskFilterIntensity2->SetNumberOfThreads( this->GetNumberOfThreads() );
        maskFilterIntensity2->SetCoordinateTolerance( m_Tol );
        maskFilterIntensity2->SetDirectionTolerance( m_Tol );
        maskFilterIntensity2->GraftOutput( this->GetOutputIntensityImage2() );
        maskFilterIntensity2->Update();
        this->GraftNthOutput( 10 , maskFilterIntensity2->GetOutput() );

        MinimumFilterTypeUC::Pointer filtermin1 = MinimumFilterTypeUC::New();
        filtermin1->SetInput1( maskFilterIntensity1->GetOutput() );
        filtermin1->SetInput2( maskFilterIntensity2->GetOutput() );
        filtermin1->SetNumberOfThreads( this->GetNumberOfThreads() );
        filtermin1->SetCoordinateTolerance( m_Tol );
        filtermin1->SetDirectionTolerance( m_Tol );

        MinimumFilterTypeUC::Pointer filtermin2 = MinimumFilterTypeUC::New();
        filtermin2->SetInput1( filtermin1->GetOutput() );
        filtermin2->SetInput2( m_LesionsDetectionImage );
        filtermin2->SetNumberOfThreads( this->GetNumberOfThreads() );
        filtermin2->SetCoordinateTolerance( m_Tol );
        filtermin2->SetDirectionTolerance( m_Tol );
        filtermin2->Update();

        OutliersIntense = filtermin2->GetOutput();
    }
    else
    {
        OutliersIntense = m_LesionsDetectionImage;
    }


    typename RemoveTouchingBorderFilterType::Pointer BorderMaskFilter = RemoveTouchingBorderFilterType::New();
    if( m_RemoveBorder )
    {
        std::cout << "Remove lesions touching border..." << std::endl;
        BorderMaskFilter->SetMask( m_MaskUC );
        BorderMaskFilter->SetInputImageSeg( OutliersIntense );
        BorderMaskFilter->SetVerbose( m_Verbose );
        BorderMaskFilter->SetTol( m_Tol );
        BorderMaskFilter->SetNumberOfThreads( this->GetNumberOfThreads() );
        BorderMaskFilter->Update();
    }

    ImageTypeUC::Pointer bigLesionsImage = ImageTypeUC::New();

    if( m_LesionSegmentationType != manualGC )
    {
        std::cout << "Computing white matter ratio..." << std::endl;
        BinaryThresholdImageFilterType_F_UC::Pointer thresholdFilterMapWM  = BinaryThresholdImageFilterType_F_UC::New();
        thresholdFilterMapWM->SetInput( m_MahalanobisFilter->GetOutputMahaWM() );
        thresholdFilterMapWM->SetLowerThreshold( m_ThresoldWMmap );
        thresholdFilterMapWM->SetInsideValue( insideValue );
        thresholdFilterMapWM->SetOutsideValue( outsideValue );
        thresholdFilterMapWM->SetNumberOfThreads( this->GetNumberOfThreads() );
        thresholdFilterMapWM->SetCoordinateTolerance( m_Tol );
        thresholdFilterMapWM->SetDirectionTolerance( m_Tol );

        CheckStructureNeighborFilterFilterType::Pointer CheckStructureNeighborFilterFilter = CheckStructureNeighborFilterFilterType::New();
        CheckStructureNeighborFilterFilter->SetInputMap( thresholdFilterMapWM->GetOutput() );
        CheckStructureNeighborFilterFilter->SetNumberOfThreads( this->GetNumberOfThreads() );
        CheckStructureNeighborFilterFilter->SetTol( m_Tol );
        if( m_RemoveBorder )
        {
            CheckStructureNeighborFilterFilter->SetInputClassification( BorderMaskFilter->GetOutputNonTouchingBorder() );
        }
        else
        {
            CheckStructureNeighborFilterFilter->SetInputClassification( OutliersIntense );
        }
        CheckStructureNeighborFilterFilter->SetRatio( m_RatioContourWM );
        CheckStructureNeighborFilterFilter->Update();
        bigLesionsImage = CheckStructureNeighborFilterFilter->GetOutput();
    }
    else
    {
        if( m_RemoveBorder )
        {
            bigLesionsImage = BorderMaskFilter->GetOutputNonTouchingBorder();
        }
        else
        {
            bigLesionsImage = OutliersIntense ;
        }
    }


    /** Size Rule: In order to avoid false positives, candidate lesions smaller than 9 mm3 in size are rejected.
     * These small candidate lesions are usually produced by noise or flow artifacts.
     * In clinical practice, lesions must have a radius of 3 mm on one image slice to be considered as such.
     */
    ImageTypeUC::SpacingType spacing = m_MaskUC->GetSpacing();
    ImageTypeUC::SpacingValueType spacingTot = spacing[0];
    for (unsigned int i = 1; i < 3;++i)
    {
        spacingTot *= spacing[i];
    }

    // Compute minimum lesion size in number of pixels
    double epsilon = 10e-6;
    double minSizeInVoxelD = m_MinLesionsSize / spacingTot;
    minSizeInVoxelD -= epsilon;
    double minSizeInVoxelD_ceil = std::ceil( minSizeInVoxelD );
    unsigned int minSizeInVoxel = static_cast<unsigned int>( minSizeInVoxelD_ceil );
    if( minSizeInVoxel > 1 )
    {
        std::cout << "Removing lesions smaller than " << minSizeInVoxel << "mm3 ..." << std::endl;
    }

    bool connectivity = false;
    ConnectedComponentType::Pointer ccFilterSize = ConnectedComponentType::New();
    ccFilterSize->SetInput( bigLesionsImage );
    ccFilterSize->SetFullyConnected( connectivity );
    ccFilterSize->SetNumberOfThreads( this->GetNumberOfThreads() );
    ccFilterSize->SetCoordinateTolerance( m_Tol );
    ccFilterSize->SetDirectionTolerance( m_Tol );

    RelabelComponentType::Pointer relabelFilterSize = RelabelComponentType::New();
    relabelFilterSize->SetInput( ccFilterSize->GetOutput() );
    relabelFilterSize->SetMinimumObjectSize( minSizeInVoxel );
    relabelFilterSize->SetNumberOfThreads( this->GetNumberOfThreads() );
    relabelFilterSize->SetCoordinateTolerance( m_Tol );
    relabelFilterSize->SetDirectionTolerance( m_Tol );
    relabelFilterSize->Update();

    m_LabeledLesions = ImageTypeInt::New();
    m_LabeledLesions->Graft(relabelFilterSize->GetOutput());
}

template <typename TInputImage>
void
GcStremMsLesionsSegmentationFilter <TInputImage>::ComputeNABT()
{
    std::cout << "Computing normal appearing brain tissues..." << std::endl;
    OutputImagePointer ImageCSF = this->GetOutputCSF();
    ImageCSF->SetRegions( m_MaskUC->GetLargestPossibleRegion() );
    ImageCSF->CopyInformation( m_MaskUC );
    ImageCSF->Allocate();
    ImageCSF->FillBuffer(0);

    OutputImagePointer ImageGM = this->GetOutputGM();
    ImageGM->SetRegions( m_MaskUC->GetLargestPossibleRegion() );
    ImageGM->CopyInformation( m_MaskUC );
    ImageGM->Allocate();
    ImageGM->FillBuffer(0);

    OutputImagePointer ImageWM = this->GetOutputWM();
    ImageWM->SetRegions( m_MaskUC->GetLargestPossibleRegion() );
    ImageWM->CopyInformation( m_MaskUC );
    ImageWM->Allocate();
    ImageWM->FillBuffer(0);

    OutputImagePointer ImageWholeSeg = this->GetOutputWholeSeg();
    ImageWholeSeg->SetRegions( m_MaskUC->GetLargestPossibleRegion() );
    ImageWholeSeg->CopyInformation( m_MaskUC );
    ImageWholeSeg->Allocate();
    ImageWholeSeg->FillBuffer(0);

    OutputImagePointer ImageLesions = this->GetOutputLesions();
    ImageLesions->SetRegions( m_MaskUC->GetLargestPossibleRegion() );
    ImageLesions->CopyInformation( m_MaskUC );
    ImageLesions->Allocate();
    ImageLesions->FillBuffer(0);

    OutputIteratorType itOutputImageCSF( ImageCSF, ImageCSF->GetLargestPossibleRegion() );
    OutputIteratorType itOutputImageGM( ImageGM, ImageGM->GetLargestPossibleRegion() );
    OutputIteratorType itOutputImageWM( ImageWM, ImageWM->GetLargestPossibleRegion() );
    OutputIteratorType itOutputImageWholeSeg( ImageWholeSeg, ImageWholeSeg->GetLargestPossibleRegion() );
    OutputIteratorType itOutputImageLesions( ImageLesions, ImageLesions->GetLargestPossibleRegion() );
    ImageIteratorTypeInt itImageLesions( m_LabeledLesions, m_LabeledLesions->GetLargestPossibleRegion() );

    ImageIteratorTypeUC itMask( m_MaskUC, m_MaskUC->GetLargestPossibleRegion() );
    while (!itMask.IsAtEnd())
    {
        if((itMask.Get()!=0) && (itImageLesions.Get() != 0))
        {
            itOutputImageLesions.Set(1);
            itOutputImageWholeSeg.Set(static_cast<OutputPixelType>(m_LabelLesions));
        }
        ++itOutputImageLesions;
        ++itOutputImageWholeSeg;
        ++itImageLesions;
        ++itMask;
    }

    itOutputImageWholeSeg.GoToBegin();
    if(m_LesionSegmentationType!=manualGC)
    {
        double minimumThWM = 0.3; // set a minimum proba value to avoid classifying grey matter as white matter

        ImageIteratorTypeD itmahaWM(m_MahalanobisFilter->GetOutputMahaWM(),m_MahalanobisFilter->GetOutputMahaWM()->GetLargestPossibleRegion());
        ImageIteratorTypeD itmahaGM(m_MahalanobisFilter->GetOutputMahaGM(),m_MahalanobisFilter->GetOutputMahaGM()->GetLargestPossibleRegion());
        ImageIteratorTypeD itmahaCSF(m_MahalanobisFilter->GetOutputMahaCSF(),m_MahalanobisFilter->GetOutputMahaCSF()->GetLargestPossibleRegion());

        itMask.GoToBegin();
        itImageLesions.GoToBegin();

        while (!itMask.IsAtEnd())
        {
            if((itMask.Get()!=0) && (itImageLesions.Get() == 0))
            {
                if( (itmahaCSF.Get()>itmahaGM.Get()) && (itmahaCSF.Get()>itmahaWM.Get()) )
                {
                    itOutputImageCSF.Set(1);
                    itOutputImageWholeSeg.Set(static_cast<OutputPixelType>(m_LabelCSF));
                }
                else if( (itmahaGM.Get()>itmahaCSF.Get()) && (itmahaGM.Get()>itmahaWM.Get()) )
                {
                    itOutputImageGM.Set(1);
                    itOutputImageWholeSeg.Set(static_cast<OutputPixelType>(m_LabelGM));
                }
                else if( (itmahaWM.Get()>itmahaGM.Get()) && (itmahaWM.Get()>itmahaCSF.Get()) && (itmahaWM.Get()>minimumThWM) )
                {
                    itOutputImageWM.Set(1);
                    itOutputImageWholeSeg.Set(static_cast<OutputPixelType>(m_LabelWM));
                }
            }
            ++itOutputImageCSF;
            ++itOutputImageGM;
            ++itOutputImageWM;
            ++itOutputImageWholeSeg;
            ++itOutputImageLesions;
            ++itImageLesions;
            ++itmahaWM;
            ++itmahaGM;
            ++itmahaCSF;
            ++itMask;
        }
    }
}


template <typename TInputImage>
void
GcStremMsLesionsSegmentationFilter <TInputImage>::GenerateData()
{
    this->CheckInputImages();

    this->RescaleImages();

    // Compute NABT model estimation
    if (m_LesionSegmentationType != manualGC)
    {
        this->ComputeAutomaticInitialization();
    }

    // Lesions detection by graph cut or thresholding
    if( m_LesionSegmentationType == strem )
    {
        this->StremThreshold();
        m_LesionsDetectionImage = this->GetOutputStrem();
    }
    else
    {
        this->GraphCut();
        m_LesionsDetectionImage = this->GetOutputGraphCut();
    }

    // Remove false positives
    this->ApplyHeuristicRules();

    // Compute NABT maps
    this->ComputeNABT();
}

}
