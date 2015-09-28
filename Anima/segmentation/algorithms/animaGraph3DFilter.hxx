#pragma once

#include "animaGraph3DFilter.h"

namespace anima
{

template <typename TInput, typename TOutput>
void Graph3DFilter<TInput, TOutput>::SetInputImage1(const TInput* image)
{
    this->SetNthInput(0, const_cast<TInput*>(image));
}

template <typename TInput, typename TOutput>
void Graph3DFilter<TInput, TOutput>::SetInputSeedProbaSources(const TSeedProba* proba)
{
    this->SetNthInput(1, const_cast<TSeedProba*>(proba));
}

template <typename TInput, typename TOutput>
void Graph3DFilter<TInput, TOutput>::SetInputSeedProbaSinks(const TSeedProba* proba)
{
    this->SetNthInput(2, const_cast<TSeedProba*>(proba));
}

template <typename TInput, typename TOutput>
void Graph3DFilter<TInput, TOutput>::SetMask(const TMask* mask)
{
    this->SetNthInput(3, const_cast<TMask*>(mask));
}

template <typename TInput, typename TOutput>
void Graph3DFilter<TInput, TOutput>::SetInputImage2(const TInput* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInput*>(image));
    m_IndexImage2 = m_NbInputs;
    m_NbInputs++;
}

template <typename TInput, typename TOutput>
void Graph3DFilter<TInput, TOutput>::SetInputImage3(const TInput* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInput*>(image));
    m_IndexImage3 = m_NbInputs;
    m_NbInputs++;
}

template <typename TInput, typename TOutput>
void Graph3DFilter<TInput, TOutput>::SetInputImage4(const TInput* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInput*>(image));
    m_IndexImage4 = m_NbInputs;
    m_NbInputs++;
}

template <typename TInput, typename TOutput>
void Graph3DFilter<TInput, TOutput>::SetInputImage5(const TInput* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInput*>(image));
    m_IndexImage5 = m_NbInputs;
    m_NbInputs++;
}


template <typename TInput, typename TOutput>
typename TInput::ConstPointer Graph3DFilter<TInput, TOutput>::GetInputImage1()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(0) );
}

template <typename TInput, typename TOutput>
itk::Image <double,3>::ConstPointer Graph3DFilter<TInput, TOutput>::GetInputSeedProbaSources()
{
    return static_cast< const TSeedProba * >
            ( this->itk::ProcessObject::GetInput(1) );
}

template <typename TInput, typename TOutput>
itk::Image <double,3>::ConstPointer Graph3DFilter<TInput, TOutput>::GetInputSeedProbaSinks()
{
    return static_cast< const TSeedProba * >
            ( this->itk::ProcessObject::GetInput(2) );
}

template <typename TInput, typename TOutput>
itk::Image <unsigned char,3>::ConstPointer Graph3DFilter<TInput, TOutput>::GetMask()
{
    return static_cast< const TMask * >
            ( this->itk::ProcessObject::GetInput(3) );
}

template <typename TInput, typename TOutput>
typename TInput::ConstPointer Graph3DFilter<TInput, TOutput>::GetInputImage2()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage2) );
}

template <typename TInput, typename TOutput>
typename TInput::ConstPointer Graph3DFilter<TInput, TOutput>::GetInputImage3()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage3) );
}

template <typename TInput, typename TOutput>
typename TInput::ConstPointer Graph3DFilter<TInput, TOutput>::GetInputImage4()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage4) );
}

template <typename TInput, typename TOutput>
typename TInput::ConstPointer Graph3DFilter<TInput, TOutput>::GetInputImage5()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage5) );
}

template <typename TInputImage, typename TOutput>
itk::DataObject::Pointer Graph3DFilter <TInputImage, TOutput>::MakeOutput(unsigned int idx)
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
typename TOutput::Pointer Graph3DFilter <TInputImage, TOutput>::GetOutput()
{
    return dynamic_cast< TOutput* >( this->itk::ProcessObject::GetOutput(0) );
}
template <typename TInputImage, typename TOutput>
typename TOutput::Pointer Graph3DFilter <TInputImage, TOutput>::GetOutputBackground()
{
    return dynamic_cast< TOutput* >( this->itk::ProcessObject::GetOutput(1) );
}


template <typename TInput, typename TOutput>
void
Graph3DFilter<TInput, TOutput>
::GenerateData()
{
    // find out if graph fit into memory
    bool mem = this->CheckMemory();

    if(mem)
    {
        // process graph with without downsampling
        this->ProcessGraphCut();
    }
    else
    {
        // if graph does not fit into memory, decimate images until we can process the graph cut.
        this->FindDownsampleFactor(); // find out how much we need to downsample

        int current_count = m_Count;

        std::cout << "-- Downsampling images for the graph... " << std::endl;
        this->InitResampleFilters();

        typename TInput::SizeType outputSize;
        typename TInput::SpacingType outputSpacing;

        while(current_count != -1)
        {
            // Resize
            outputSize[0] = m_InputSize[0] / m_DownsamplingFactor;
            outputSize[1] = m_InputSize[1] / m_DownsamplingFactor;
            outputSize[2] = m_InputSize[2] / m_DownsamplingFactor;

            outputSpacing[0] = this->GetMask()->GetSpacing()[0] * (static_cast<double>(m_InputSize[0]) / static_cast<double>(outputSize[0]));
            outputSpacing[1] = this->GetMask()->GetSpacing()[1] * (static_cast<double>(m_InputSize[1]) / static_cast<double>(outputSize[1]));
            outputSpacing[2] = this->GetMask()->GetSpacing()[2] * (static_cast<double>(m_InputSize[2]) / static_cast<double>(outputSize[2]));

            m_ResampleMask->SetSize(outputSize);
            m_ResampleMask->SetOutputSpacing(outputSpacing);
            m_ResampleMask->Update();

            m_ResampleSources->SetSize(outputSize);
            m_ResampleSources->SetOutputSpacing(outputSpacing);

            m_ResampleSinks->SetSize(outputSize);
            m_ResampleSinks->SetOutputSpacing(outputSpacing);

            m_Resample1->SetSize(outputSize);
            m_Resample1->SetOutputSpacing(outputSpacing);

            m_Resample2->SetSize(outputSize);
            m_Resample2->SetOutputSpacing(outputSpacing);

            m_Resample3->SetSize(outputSize);
            m_Resample3->SetOutputSpacing(outputSpacing);

            m_Resample4->SetSize(outputSize);
            m_Resample4->SetOutputSpacing(outputSpacing);

            m_Resample5->SetSize(outputSize);
            m_Resample5->SetOutputSpacing(outputSpacing);

            if(current_count == m_Count)
            {
                m_CurrentMask = m_ResampleMask->GetOutput();
            }

            // Process graph cut with downsampled images
            this->ProcessDownsampledGraphCut(current_count);

            // Process the graph cut with the following downsampling factor
            m_DownsamplingFactor/=2;

            typename TInput::SizeType upSize;
            upSize[0] = m_InputSize[0] / m_DownsamplingFactor;
            upSize[1] = m_InputSize[1] / m_DownsamplingFactor;
            upSize[2] = m_InputSize[2] / m_DownsamplingFactor;

            if(current_count > 0)
            {
                TMask::Pointer upImage = Upsample(m_NLinksFilterDecim->GetOutput(), upSize, m_OutputDirection, m_OutputOrigin);

                typename TMask::SpacingType maskSpacing;
                maskSpacing[0] = upImage->GetSpacing()[0];
                maskSpacing[1] = upImage->GetSpacing()[1];
                maskSpacing[2] = upImage->GetSpacing()[2];

                m_ResampleMask->SetSize(upSize);
                m_ResampleMask->SetOutputSpacing(maskSpacing);
                m_ResampleMask->Update();

                // Dilate image inside the mask
                m_CurrentMask = Dilate(upImage, m_ResampleMask->GetOutput());
            }

            current_count--;
        }
    }
}

template <typename TInput, typename TOutput>
bool
Graph3DFilter<TInput, TOutput>
::CheckMemory()
{
    unsigned int nb_vox = 0;
    bool mem = true;
    MaskRegionConstIteratorType maskIt (this->GetMask(),this->GetMask()->GetLargestPossibleRegion() );
    while (!maskIt.IsAtEnd())
    {
        if (maskIt.Get() != 0)
        {
            nb_vox++;
        }
        ++maskIt;
    }

    unsigned int nb_edges = 7*nb_vox;

    try
    {
        m_graph = new GraphType(nb_vox, nb_edges);
    }
    catch (std::bad_alloc& ba)
    {
        std::cerr << "-- In Graph3DFilter: insufficient memory to create the graph: " << ba.what() << '\n';
        mem = false;
    }

    if (m_graph)
    {
        delete m_graph;
        m_graph = NULL;
    }
    return mem;
}

template <typename TInput, typename TOutput>
void
Graph3DFilter<TInput, TOutput>
::ProcessGraphCut()
{
    m_NLinksFilter->SetSigma( this->GetSigma() );
    m_NLinksFilter->SetUseSpectralGradient( this->GetUseSpectralGradient() );
    m_NLinksFilter->SetMatrix( this->GetMatrix() );
    m_NLinksFilter->SetMatFilename( this->GetMatFilename() );
    m_NLinksFilter->SetVerbose( this->GetVerbose() );
    m_NLinksFilter->SetNumberOfThreads(this->GetNumberOfThreads());
    m_NLinksFilter->SetTol( m_Tol );

    m_NLinksFilter->SetInputSeedProbaSources( this->GetInputSeedProbaSources() );
    m_NLinksFilter->SetInputSeedProbaSinks( this->GetInputSeedProbaSinks() );
    m_NLinksFilter->SetMask( this->GetMask() );

    m_NLinksFilter->SetInputImage1( this->GetInputImage1() );
    if(m_IndexImage2!=m_NbMaxImages){m_NLinksFilter->SetInputImage2( this->GetInputImage2() );}
    if(m_IndexImage3!=m_NbMaxImages){m_NLinksFilter->SetInputImage3( this->GetInputImage3() );}
    if(m_IndexImage4!=m_NbMaxImages){m_NLinksFilter->SetInputImage4( this->GetInputImage4() );}
    if(m_IndexImage5!=m_NbMaxImages){m_NLinksFilter->SetInputImage5( this->GetInputImage5() );}

    m_NLinksFilter->GraftNthOutput( 0, this->GetOutput() );
    m_NLinksFilter->GraftNthOutput( 1, this->GetOutputBackground() );

    try
    {
        m_NLinksFilter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        exit(-1);
    }

    this->GraftNthOutput( 0 , m_NLinksFilter->GetOutput() );
    this->GraftNthOutput( 1 , m_NLinksFilter->GetOutputBackground() );
}

template <typename TInput, typename TOutput>
void
Graph3DFilter<TInput, TOutput>::FindDownsampleFactor()
{
    m_DownsamplingFactor = 2.0;
    m_Count = 1;

    bool mem2 = false;
    while(!mem2)
    {
        typename TInput::SizeType inputSize = this->GetMask()->GetLargestPossibleRegion().GetSize();
        typename TInput::DirectionType outputDirection = this->GetMask()->GetDirection();
        typename TInput::PointType outputOrigin = this->GetMask()->GetOrigin();

        // Resize images
        typename TInput::SizeType outputSize;
        outputSize[0] = inputSize[0] / m_DownsamplingFactor;
        outputSize[1] = inputSize[1] / m_DownsamplingFactor;
        outputSize[2] = inputSize[2] / m_DownsamplingFactor;
        typename TInput::SpacingType outputSpacing;
        outputSpacing[0] = this->GetMask()->GetSpacing()[0] * (static_cast<double>(inputSize[0]) / static_cast<double>(outputSize[0]));
        outputSpacing[1] = this->GetMask()->GetSpacing()[1] * (static_cast<double>(inputSize[1]) / static_cast<double>(outputSize[1]));
        outputSpacing[2] = this->GetMask()->GetSpacing()[2] * (static_cast<double>(inputSize[2]) / static_cast<double>(outputSize[2]));

        ResampleImageFilterMaskType::Pointer resampleMask = ResampleImageFilterMaskType::New();
        resampleMask->SetOutputDirection( outputDirection );
        resampleMask->SetOutputOrigin( outputOrigin );
        resampleMask->SetInput( this->GetMask() );
        resampleMask->SetSize( outputSize );
        resampleMask->SetOutputSpacing( outputSpacing );
        resampleMask->SetTransform( TransformType::New() );
        resampleMask->SetNumberOfThreads( this->GetNumberOfThreads() );
        resampleMask->SetCoordinateTolerance( m_Tol );
        resampleMask->SetDirectionTolerance( m_Tol );
        resampleMask->Update();

        unsigned int nb_vox = 0;
        MaskRegionIteratorType maskIt (resampleMask->GetOutput(),resampleMask->GetOutput()->GetLargestPossibleRegion() );
        while (!maskIt.IsAtEnd())
        {
            if (maskIt.Get() != 0)
            {
                nb_vox++;
            }
            ++maskIt;
        }

        unsigned int nb_edges = 7*nb_vox;
        mem2 = true;
        try
        {
            m_graph = new GraphType(nb_vox, nb_edges);
        }
        catch (std::bad_alloc& ba)
        {
            std::cerr << "-- In Graph3DFilter: insufficient memory to create the graph ICI: " << ba.what() << '\n';
            mem2 = false;
            m_Count++;
            m_DownsamplingFactor*=2.0;
        }

        if (m_graph)
        {
            delete m_graph;
            m_graph = NULL;
        }
    }
}

template <typename TInput, typename TOutput>
void
Graph3DFilter<TInput, TOutput>
::InitResampleFilters()
{
    m_InputSize = this->GetMask()->GetLargestPossibleRegion().GetSize();
    m_OutputDirection = this->GetMask()->GetDirection();
    m_OutputOrigin = this->GetMask()->GetOrigin();

    m_ResampleMask = ResampleImageFilterMaskType::New();
    m_ResampleMask->SetOutputDirection(m_OutputDirection);
    m_ResampleMask->SetOutputOrigin(m_OutputOrigin);
    m_ResampleMask->SetInput(this->GetMask());
    m_ResampleMask->SetTransform(TransformType::New());
    m_ResampleMask->SetNumberOfThreads(this->GetNumberOfThreads());
    m_ResampleMask->SetCoordinateTolerance(m_Tol);
    m_ResampleMask->SetDirectionTolerance(m_Tol);

    m_ResampleSources = ResampleImageFilterFloatType::New();
    m_ResampleSources->SetOutputDirection(m_OutputDirection);
    m_ResampleSources->SetOutputOrigin(m_OutputOrigin);
    m_ResampleSources->SetInput(this->GetInputSeedProbaSources());
    m_ResampleSources->SetTransform(TransformType::New());
    m_ResampleSources->SetNumberOfThreads(this->GetNumberOfThreads());
    m_ResampleSources->SetCoordinateTolerance(m_Tol);
    m_ResampleSources->SetDirectionTolerance(m_Tol);

    m_ResampleSinks = ResampleImageFilterFloatType::New();
    m_ResampleSinks->SetOutputDirection(m_OutputDirection);
    m_ResampleSinks->SetOutputOrigin(m_OutputOrigin);
    m_ResampleSinks->SetInput(this->GetInputSeedProbaSinks());
    m_ResampleSinks->SetTransform(TransformType::New());
    m_ResampleSinks->SetNumberOfThreads(this->GetNumberOfThreads());
    m_ResampleSinks->SetCoordinateTolerance(m_Tol);
    m_ResampleSinks->SetDirectionTolerance(m_Tol);

    m_Resample1 = ResampleImageFilterType::New();
    m_Resample1->SetOutputDirection(m_OutputDirection);
    m_Resample1->SetOutputOrigin(m_OutputOrigin);
    m_Resample1->SetInput(this->GetInputImage1());
    m_Resample1->SetTransform(TransformType::New());
    m_Resample1->SetNumberOfThreads(this->GetNumberOfThreads());
    m_Resample1->SetCoordinateTolerance(m_Tol);
    m_Resample1->SetDirectionTolerance(m_Tol);

    m_Resample2 = ResampleImageFilterType::New();
    m_Resample2->SetOutputDirection(m_OutputDirection);
    m_Resample2->SetOutputOrigin(m_OutputOrigin);
    m_Resample2->SetInput(this->GetInputImage2());
    m_Resample2->SetTransform(TransformType::New());
    m_Resample2->SetNumberOfThreads(this->GetNumberOfThreads());
    m_Resample2->SetCoordinateTolerance(m_Tol);
    m_Resample2->SetDirectionTolerance(m_Tol);

    m_Resample3 = ResampleImageFilterType::New();
    m_Resample3->SetOutputDirection(m_OutputDirection);
    m_Resample3->SetOutputOrigin(m_OutputOrigin);
    m_Resample3->SetInput(this->GetInputImage4());
    m_Resample3->SetTransform(TransformType::New());
    m_Resample3->SetNumberOfThreads(this->GetNumberOfThreads());
    m_Resample3->SetCoordinateTolerance(m_Tol);
    m_Resample3->SetDirectionTolerance(m_Tol);

    m_Resample4 = ResampleImageFilterType::New();
    m_Resample4->SetOutputDirection(m_OutputDirection);
    m_Resample4->SetOutputOrigin(m_OutputOrigin);
    m_Resample4->SetInput(this->GetInputImage4());
    m_Resample4->SetTransform(TransformType::New());
    m_Resample4->SetNumberOfThreads(this->GetNumberOfThreads());
    m_Resample4->SetCoordinateTolerance(m_Tol);
    m_Resample4->SetDirectionTolerance(m_Tol);

    m_Resample5 = ResampleImageFilterType::New();
    m_Resample5->SetOutputDirection(m_OutputDirection);
    m_Resample5->SetOutputOrigin(m_OutputOrigin);
    m_Resample5->SetInput(this->GetInputImage5());
    m_Resample5->SetTransform(TransformType::New());
    m_Resample5->SetNumberOfThreads(this->GetNumberOfThreads());
    m_Resample5->SetCoordinateTolerance(m_Tol);
    m_Resample5->SetDirectionTolerance(m_Tol);
}

template <typename TInput, typename TOutput>
void
Graph3DFilter<TInput, TOutput>
::ProcessDownsampledGraphCut(int current_count)
{
    m_NLinksFilterDecim = NLinksFilterType::New();
    m_NLinksFilterDecim->SetSigma( this->GetSigma() );
    m_NLinksFilterDecim->SetUseSpectralGradient( this->GetUseSpectralGradient() );
    m_NLinksFilterDecim->SetMatrix( this->GetMatrix() );
    m_NLinksFilterDecim->SetMatFilename( this->GetMatFilename() );
    m_NLinksFilterDecim->SetVerbose( this->GetVerbose() );
    m_NLinksFilterDecim->SetNumberOfThreads(this->GetNumberOfThreads());
    m_NLinksFilterDecim->SetTol( m_Tol );
    m_NLinksFilterDecim->SetMask( m_CurrentMask ); // mandatory brain mask

    m_NLinksFilterDecim->SetInputSeedProbaSources( m_ResampleSources->GetOutput() );
    m_NLinksFilterDecim->SetInputSeedProbaSinks( m_ResampleSinks->GetOutput() );

    m_NLinksFilterDecim->SetInputImage1( m_Resample1->GetOutput() );
    if(m_IndexImage2!=m_NbMaxImages){m_NLinksFilterDecim->SetInputImage2( m_Resample2->GetOutput() );}
    if(m_IndexImage3!=m_NbMaxImages){m_NLinksFilterDecim->SetInputImage3( m_Resample3->GetOutput() );}
    if(m_IndexImage4!=m_NbMaxImages){m_NLinksFilterDecim->SetInputImage4( m_Resample4->GetOutput() );}
    if(m_IndexImage5!=m_NbMaxImages){m_NLinksFilterDecim->SetInputImage5( m_Resample5->GetOutput() );}

    if(current_count==0)
    {
        m_NLinksFilterDecim->GraftNthOutput( 0, this->GetOutput() );
        m_NLinksFilterDecim->GraftNthOutput( 1, this->GetOutputBackground() );
    }

    try
    {
        m_NLinksFilterDecim->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        exit(-1);
    }

    if(current_count==0)
    {
        this->GraftNthOutput( 0 , m_NLinksFilterDecim->GetOutput() );
        this->GraftNthOutput( 1 , m_NLinksFilterDecim->GetOutputBackground() );
    }
}

template <typename TInput, typename TOutput>
itk::Image <unsigned char,3>::Pointer
Graph3DFilter<TInput, TOutput>
::Dilate(const TMask* input, const TMask* mask)
{
    unsigned int radius = 4;

    typedef itk::BinaryBallStructuringElement<TMask::PixelType, TMask::ImageDimension> StructuringElementType;
    StructuringElementType structuringElement;
    structuringElement.SetRadius(radius);
    structuringElement.CreateStructuringElement();

    typedef itk::GrayscaleDilateImageFilter <TMask,TMask,StructuringElementType> DilateFilterType;
    DilateFilterType::Pointer dilateFilter = DilateFilterType::New();
    dilateFilter->SetInput(input);
    dilateFilter->SetKernel(structuringElement);
    dilateFilter->SetNumberOfThreads(this->GetNumberOfThreads());
    dilateFilter->SetCoordinateTolerance(m_Tol);
    dilateFilter->SetDirectionTolerance(m_Tol);

    typedef itk::MaskImageFilter< TMask, TMask > MaskFilterType;
    MaskFilterType::Pointer maskFilter = MaskFilterType::New();
    maskFilter->SetInput(dilateFilter->GetOutput());
    maskFilter->SetMaskImage(mask);
    maskFilter->Update();

    return maskFilter->GetOutput();
}

template <typename TInput, typename TOutput>
itk::Image <unsigned char,3>::Pointer
Graph3DFilter<TInput, TOutput>
::Upsample(const TMask* input, TMask::SizeType upSize, TMask::DirectionType outputDirection, TMask::PointType outputOrigin)
{
    // We don't want any transform on our image except rescaling which is not
    // specified by a transform but by the input/output spacing.
    // So no transform will be specified.
    typedef itk::IdentityTransform<double, 3> T_Transform;
    // Instantiate the transform and specify it should be the id transform.
    T_Transform::Pointer Transform = T_Transform::New();
    Transform->SetIdentity();

    // If ITK resampler determines there is something to interpolate which is
    // usually the case when upscaling (!) then we must specify the interpolation
    // algorithm. In our case, we want bicubic interpolation. One way to implement
    // it is with a third order b-spline. So the type is specified here and the
    // order will be specified with a method call later on.

    typedef itk::BSplineInterpolateImageFunction<TMask, double, double> T_Interpolator;
    // Instantiate the b-spline interpolator and set it as the third orderfor bicubic.
    T_Interpolator::Pointer Interpolator = T_Interpolator::New();
    Interpolator->SetSplineOrder(3);

    // The resampler type itself.
    typedef itk::ResampleImageFilter<TMask, TMask> T_ResampleFilter;
    // Instantiate the resampler. Wire in the transform and the interpolator.
    T_ResampleFilter::Pointer ResizeFilter = T_ResampleFilter::New();
    ResizeFilter->SetTransform(Transform);
    ResizeFilter->SetInterpolator(Interpolator);

    // Set the output origin.
    ResizeFilter->SetOutputOrigin(m_OutputOrigin);
    ResizeFilter->SetOutputDirection(m_OutputDirection);

    //     Compute and set the output spacing
    //     Compute the output spacing from input spacing and old and new sizes.
    //
    //     The computation must be so that the following holds:
    //
    //     new width         old x spacing
    //     ----------   =   ---------------
    //     old width         new x spacing
    //
    //     we specify new height and width and compute new spacings

    unsigned int nNewSize_0 = upSize[0];
    unsigned int nNewSize_1 = upSize[1];
    unsigned int nNewSize_2 = upSize[2];

    // Fetch original image size.
    const TMask::RegionType& inputRegion = input->GetLargestPossibleRegion();
    const TMask::SizeType& vnInputSize = inputRegion.GetSize();
    unsigned int nOldSize_0 = vnInputSize[0];
    unsigned int nOldSize_1 = vnInputSize[1];
    unsigned int nOldSize_2 = vnInputSize[2];

    // Fetch original image spacing.
    const TMask::SpacingType& vfInputSpacing = input->GetSpacing();

    double vfOutputSpacing[3];
    vfOutputSpacing[0] = vfInputSpacing[0] * (double) nOldSize_0 / (double) nNewSize_0;
    vfOutputSpacing[1] = vfInputSpacing[1] * (double) nOldSize_1 / (double) nNewSize_1;
    vfOutputSpacing[2] = vfInputSpacing[2] * (double) nOldSize_2 / (double) nNewSize_2;

    // Set the output spacing.
    ResizeFilter->SetOutputSpacing(vfOutputSpacing);

    // Set the output size
    itk::Size<3> vnOutputSize = { {nNewSize_0, nNewSize_1, nNewSize_2} };
    ResizeFilter->SetSize(vnOutputSize);
    ResizeFilter->SetInput(input);
    ResizeFilter->SetCoordinateTolerance(m_Tol);
    ResizeFilter->SetDirectionTolerance(m_Tol);
    ResizeFilter->SetNumberOfThreads(this->GetNumberOfThreads());
    ResizeFilter->Update();

    return ResizeFilter->GetOutput();
}

} //end of namespace anima
