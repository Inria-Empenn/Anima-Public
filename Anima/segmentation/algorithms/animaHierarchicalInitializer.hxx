#include "animaHierarchicalInitializer.h"

namespace anima
{

template <typename TInputImage, typename TMaskImage>
void HierarchicalInitializer<TInputImage,TMaskImage>::SetMask(const TMaskImage* mask)
{
    this->SetNthInput(0, const_cast<TMaskImage*>(mask));
}

template <typename TInputImage, typename TMaskImage>
void HierarchicalInitializer<TInputImage,TMaskImage>::SetInputImage1(const TInputImage* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInputImage*>(image));
    m_IndexImage1 = m_NbInputs;
    m_NbInputs++;
}

template <typename TInputImage, typename TMaskImage>
void HierarchicalInitializer<TInputImage,TMaskImage>::SetInputImage2(const TInputImage* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInputImage*>(image));
    m_IndexImage2 = m_NbInputs;
    m_NbInputs++;
}

template <typename TInputImage, typename TMaskImage>
void HierarchicalInitializer<TInputImage,TMaskImage>::SetInputImage3(const TInputImage* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInputImage*>(image));
    m_IndexImage3 = m_NbInputs;
    m_NbInputs++;
}
template <typename TInputImage, typename TMaskImage>
void HierarchicalInitializer<TInputImage,TMaskImage>::SetInputImage4(const TInputImage* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInputImage*>(image));
    m_IndexImage4 = m_NbInputs;
    m_NbInputs++;
}

template <typename TInputImage, typename TMaskImage>
void HierarchicalInitializer<TInputImage,TMaskImage>::SetInputImage5(const TInputImage* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInputImage*>(image));
    m_IndexImage5 = m_NbInputs;
    m_NbInputs++;
}


template <typename TInputImage, typename TMaskImage>
typename TMaskImage::ConstPointer HierarchicalInitializer<TInputImage,TMaskImage>::GetMask()
{
    return static_cast< const TMaskImage * >
            ( this->itk::ProcessObject::GetInput(0) );
}

template <typename TInputImage, typename TMaskImage>
typename TInputImage::ConstPointer HierarchicalInitializer<TInputImage,TMaskImage>::GetInputImage1()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage1) );
}

template <typename TInputImage, typename TMaskImage>
typename TInputImage::ConstPointer HierarchicalInitializer<TInputImage,TMaskImage>::GetInputImage2()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage2) );
}

template <typename TInputImage, typename TMaskImage>
typename TInputImage::ConstPointer HierarchicalInitializer<TInputImage,TMaskImage>::GetInputImage3()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage3) );
}

template <typename TInputImage, typename TMaskImage>
typename TInputImage::ConstPointer HierarchicalInitializer<TInputImage,TMaskImage>::GetInputImage4()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage4) );
}

template <typename TInputImage, typename TMaskImage>
typename TInputImage::ConstPointer HierarchicalInitializer<TInputImage,TMaskImage>::GetInputImage5()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage5) );
}

template <typename TInputImage, typename TMaskImage>
void HierarchicalInitializer<TInputImage,TMaskImage>
::ComputeSolution1D()
{
    // Compute EM for T1-w
    std::vector<double> min1D( 1, 255.0 ); //min values for random Initialization
    std::vector<double> max1D( 1, 0.0 ); //max values for random Initialization

    // Decide initializer
    RandomInitializer::Pointer initiaRandom = RandomInitializer::New();
    initiaRandom->SetMinValues( min1D );
    initiaRandom->SetMaxValues( max1D );
    initiaRandom->SetDimensionGaussian(1);
    initiaRandom->SetNbGaussian(3);

    typename GaussianREMEstimatorType::Pointer estimator= GaussianREMEstimatorType::New();
    estimator ->SetMaxIterations( 100 );
    estimator ->SetModelMinDistance( 1e-3 );
    estimator ->SetMaxIterationsConc( 10 );
    estimator ->SetStremMode( false );
    estimator ->SetRejectionRatio( m_Robust );
    estimator ->SetMask( this->GetMask() );
    estimator ->SetInputImage1( this->GetInputImage1() );
    estimator ->SetVerbose( false );

    std::vector<unsigned int> emSteps( 1, 60 );
    std::vector<unsigned int> iterSteps( 1, 60 );

    typename ClassificationStrategyType::Pointer strategy = ClassificationStrategyType::New();
    strategy->SetEstimator( estimator );
    strategy->SetInitializer( initiaRandom );
    strategy->SetStrategy( emSteps, iterSteps );
    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback ->SetCallback(eventCallback);
    strategy ->AddObserver(itk::ProgressEvent(), callback );
    strategy->Update();
    std::cout << std::endl;

    ModelMap solutionClassifier;
    ModelMap2 solutionClassifierAlpha;

    if ( !strategy->GetSolutionMap( solutionClassifier, solutionClassifierAlpha ) )
    {
        std::cerr << "-- Error: No solution found in first modality" << std::endl;
        return;
    }

    m_Solution1D = solutionClassifier.begin() ->second;
    m_Solution1DAlphas = solutionClassifierAlpha.begin() ->second;
}

template <typename TInputImage, typename TMaskImage>
void HierarchicalInitializer<TInputImage,TMaskImage>
::ComputeInitialT1Classification( )
{
    // Classify T1-w
    typename ImageClassifierType::Pointer classifier = ImageClassifierType::New();
    classifier->SetInputImage1(this->GetInputImage1());
    classifier->SetMask(this->GetMask());
    classifier->SetGaussianModel(m_Solution1D);
    classifier->SetAlphas(m_Solution1DAlphas);
    classifier->SetTol(m_Tol);
    classifier->Update();

    // creating one image per class (CSF, GM, WM)
    std::vector<InputIteratorType> classesIt;
    for(int i = 0; i < 3; i++)
    {
        InputImagePointer im= InputImageType::New();
        im->SetRegions(this->GetMask()->GetLargestPossibleRegion());
        im->CopyInformation(this->GetMask());
        im->Allocate();
        im->FillBuffer(0);
        m_ImagesClasses.push_back(im);

        InputIteratorType imIt(im, im->GetLargestPossibleRegion() );
        classesIt.push_back(imIt);
    }

    InputConstIteratorType classifIt(classifier->GetOutput(), classifier->GetOutput()->GetLargestPossibleRegion() );
    while(!classifIt.IsAtEnd())
    {
        if(!(classifIt.Get() == 0 || classifIt.Get() > 3 ))
        {
            classesIt[ classifIt.Get() - 1 ].Set(1);
        }
        ++classifIt;
        for(unsigned int i = 0; i < classesIt.size(); i++)
        {
            ++classesIt[i];
        }
    }
}

template <typename TInputImage, typename TMaskImage>
void HierarchicalInitializer<TInputImage,TMaskImage>
::ComputeHistogram( )
{
    const double resolution = 1;
    const double bandwidth = 5;

    // selection rules [sequence][class]
    // +1 brighter max
    // 0 absolute max
    // -1 darker max
    std::vector< std::vector<int> > selectionRules( 3, std::vector<int>( 3, 0 ) );
    selectionRules[ 1 ][ 0 ] = 1; // T2-w & CSF: brighter max
    if(m_ThirdIsFLAIR)
        selectionRules[ 2 ][ 0 ] = -1; // FLAIR & CSF: darker max
    else
        selectionRules[ 2 ][ 0 ] = 1; // PD-w & CSF: brighter max

    m_SelectedMax.resize( 3, std::vector<double>( 3, 0.0 ) );
    m_Stds.resize(3, std::vector<double>( 3, 0.0 ));

    typedef double MeasurementType;
    const unsigned int MeasurementVectorLength = 1;
    typedef itk::Vector< MeasurementType , MeasurementVectorLength > MeasurementVectorType;
    typedef itk::Statistics::ListSample< MeasurementVectorType > ListSampleType;
    typedef itk::Statistics::DenseFrequencyContainer2 FrequencyContainerType;
    typedef itk::Statistics::Histogram< MeasurementType,FrequencyContainerType > HistogramType;
    typedef itk::Statistics::SampleToHistogramFilter<ListSampleType,HistogramType> FilterType;

    InputImagePixelType valMaxInput = std::numeric_limits<InputImagePixelType>::max();
    InputImagePixelType valMinInput = std::numeric_limits<InputImagePixelType>::min();

    // Compute mean values
    for ( unsigned int im = 1; im < m_NumberOfModalities; im++ )
    {
        for ( unsigned int c = 0; c < m_ImagesClasses.size(); c++ )
        {
            //create the histogram with the selected mask
            typename MaskFilterType::Pointer maskFilter = MaskFilterType::New();
            typename MinMaxCalculatorType::Pointer calculator = MinMaxCalculatorType::New();
            maskFilter->SetDirectionTolerance(m_Tol);
            maskFilter->SetCoordinateTolerance(m_Tol);
            maskFilter->SetInput( m_ImagesVector[ im ] );
            maskFilter->SetMaskImage( m_ImagesClasses[ c ] );

            maskFilter->SetOutsideValue( valMaxInput );
            maskFilter->Update();
            calculator->SetImage( maskFilter->GetOutput() );
            calculator->ComputeMinimum();
            InputImagePixelType minimumResult = calculator->GetMinimum();

            maskFilter->SetOutsideValue( valMinInput );
            maskFilter->Update();
            calculator->SetImage( maskFilter->GetOutput() );
            calculator->ComputeMaximum();
            InputImagePixelType maximumResult = calculator->GetMaximum();

            double minV = static_cast<double>(minimumResult);
            double maxV = static_cast<double>(maximumResult);

            unsigned int bins = static_cast<unsigned int>(std::floor( maxV - minV + 1.0 ));

            ListSampleType::Pointer listSample = ListSampleType::New();
            listSample->SetMeasurementVectorSize( MeasurementVectorLength );

            InputConstIteratorType tmpIt (m_ImagesVector[im], m_ImagesVector[im]->GetLargestPossibleRegion() );
            InputIteratorType classesCIt (m_ImagesClasses[ c ], m_ImagesClasses[ c ]->GetLargestPossibleRegion() );

            MeasurementVectorType mv;
            while (!classesCIt.IsAtEnd())
            {
                if(classesCIt.Get()!=0)
                {
                    mv[0] = ( MeasurementType ) tmpIt.Get();
                    listSample->PushBack(mv);
                }
                ++classesCIt;
                ++tmpIt;
            }

            // creation histo
            const unsigned int numberOfComponents = 1;
            HistogramType::SizeType size( numberOfComponents );
            size.Fill( bins );
            HistogramType::MeasurementVectorType min( numberOfComponents );
            HistogramType::MeasurementVectorType max( numberOfComponents );
            min.Fill( minV );
            max.Fill( maxV );

            FilterType::Pointer filter = FilterType::New();
            filter->SetInput( listSample );
            filter->SetHistogramSize( size );
            filter->SetHistogramBinMinimum( min );
            filter->SetHistogramBinMaximum( max );
            filter->Update();

            HistogramType::ConstPointer histogram = filter->GetOutput();

            double bins2 = static_cast<double>(bins) / resolution;
            unsigned int smoothedBins = static_cast<unsigned int>(bins2 + 1);

            HistogramType::SizeType size2(MeasurementVectorLength);
            size2.Fill(smoothedBins);

            HistogramType::MeasurementVectorType lowerBound;
            lowerBound.SetSize(smoothedBins);
            lowerBound.Fill(minV);

            HistogramType::MeasurementVectorType upperBound;
            upperBound.SetSize(smoothedBins);
            upperBound.Fill(maxV);

            HistogramType::Pointer smoothed = HistogramType::New();
            smoothed->SetMeasurementVectorSize(MeasurementVectorLength);
            smoothed->Initialize(size2, lowerBound, upperBound );
            smoothed->SetFrequency(0);

            std::vector<double> smoothedFrequency(smoothedBins,0);

            for (unsigned int i = 0; i < bins; ++i)
            {
                for (unsigned int j = 0; j < smoothedBins; ++j)
                {
                    // Calculate distance in intensity values
                    double valueI = static_cast<double>(i) / static_cast<double>(bins -1) * (maxV - minV) + minV;
                    double valueJ = static_cast<double>(j) / static_cast<double>(smoothedBins -1) * (maxV - minV) + minV;
                    double kvalue = (valueI - valueJ) / bandwidth;
                    double norm = std::exp( -0.5 * kvalue * kvalue ) / std::sqrt( 2.0 * M_PI );
                    double value = norm * static_cast<double>(histogram->GetFrequency(i));
                    smoothedFrequency[j] += value;
                }
            }

            for (unsigned int i = 0; i < smoothedBins; i++ )
            {
                smoothedFrequency[i]/=(static_cast<double>( histogram->GetTotalFrequency()) * bandwidth);
            }

            // get maxs
            std::vector<unsigned int> maxs;
            const unsigned int localMax = 3;
            for (unsigned int i = 0; i < smoothedBins; i++ )
            {
		int testVal = i - localMax;
                unsigned int jmin = testVal < 0 ? 0 : testVal;	

                unsigned int jmax = ( i + localMax ) > smoothedBins ? smoothedBins : ( i + localMax );

                bool ismax = true;
                for (unsigned int j = jmin; j < jmax; j++ )
                {
                    if ( (i != j) && (smoothedFrequency[i] <= smoothedFrequency[j]) )
                    {
                        ismax = false;
                        continue;
                    }
                }

                if ( ismax && ( smoothedFrequency[i] > 0.001 ) )
                {
                    maxs.push_back( i );
                }
            }

            if(maxs.size()==0)
            {
                std::cerr << "--Error in Hierarchical Initializer: no maximum value" << std::endl;
                exit(-1);
            }
            else
            {
                // Keeping only the good max according to the rules
                double value = 0;
                switch ( selectionRules[ im ][ c ] )
                {
                case - 1: // darker max
                {
                    value = ( static_cast<double>(maxs[0]) / static_cast<double>(bins-1) ) * (maxV - minV) + minV;
                    m_SelectedMax[im][c] = value;
                    break;
                }
                case 1: // brigther max
                {
                    value = ( static_cast<double>(maxs[ maxs.size() - 1 ]) / static_cast<double>(bins-1) ) * (maxV - minV) + minV;
                    m_SelectedMax[ im ][ c ] = value;
                    break;
                }
                case 0: // absolute max
                default:
                {
                    unsigned int absMaxIndex = 0;
                    unsigned int absMaxValue = 0;
                    for ( unsigned int i = 0; i < maxs.size(); i++ )
                    {
                        if ( maxs[ i ] > absMaxValue )
                        {
                            absMaxValue = maxs[ i ];
                            absMaxIndex = i;
                        }
                    }
                    value = ( static_cast<double>(maxs[ absMaxIndex ]) / static_cast<double>(bins-1) ) * (maxV - minV) + minV;
                    m_SelectedMax[ im ][ c ] = value;
                    break;
                }
                }
            }
        }
    }
}

template <typename TInputImage, typename TMaskImage>
void HierarchicalInitializer<TInputImage,TMaskImage>
::ComputeVariances()
{
    ImageTypeF::Pointer diff = ImageTypeF::New();
    diff->SetRegions(this->GetMask()->GetLargestPossibleRegion());
    diff->CopyInformation(this->GetMask());
    diff->Allocate();
    diff->FillBuffer(0);
    ImageIteratorTypeF diffIt(diff, diff->GetLargestPossibleRegion() );

    std::vector<InputIteratorType> classesItVec;
    for(unsigned int i = 0; i < m_ImagesClasses.size(); i++)
    {
        InputIteratorType classesIt (m_ImagesClasses[i],m_ImagesClasses[i]->GetLargestPossibleRegion() );
        classesItVec.push_back(classesIt);
    }

    std::vector<InputConstIteratorType> imagesVectorItVec;
    for(unsigned int i = 0; i < m_NumberOfModalities; i++)
    {
        InputConstIteratorType imagesVectorIt (m_ImagesVector[i], m_ImagesVector[i]->GetLargestPossibleRegion() );
        imagesVectorItVec.push_back(imagesVectorIt);
    }

    const float estimatorFactor = 1.4918; // factor to compute standard deviation using a robust variance estimator

    for ( unsigned int im = 1; im < m_NumberOfModalities; im++ )
    {
        for ( unsigned int c = 0; c < m_ImagesClasses.size(); c++ )
        {

            for(unsigned int i = 0; i < m_ImagesClasses.size(); i++)
            {
                classesItVec[i].GoToBegin();
            }
            for(unsigned int i = 0; i < m_NumberOfModalities; i++)
            {
                imagesVectorItVec[i].GoToBegin();
            }
            diffIt.GoToBegin();

            while (!diffIt.IsAtEnd())
            {

                if(classesItVec[c].Get() > 0)
                {
                    diffIt.Set(std::fabs(imagesVectorItVec[im].Get()-m_SelectedMax[im][c]));
                }

                for(unsigned int i = 0; i < m_ImagesClasses.size(); i++)
                {
                    ++classesItVec[i];
                }
                for(unsigned int i = 0; i < m_NumberOfModalities; i++)
                {
                    ++imagesVectorItVec[i];
                }
                ++diffIt;
            }
            m_Stds[im][c] = estimatorFactor * regionMedianValue(diff,m_ImagesClasses[c]);
        }
    }
}

template <typename TInputImage, typename TMaskImage>
void HierarchicalInitializer<TInputImage,TMaskImage>
::FillInitialGaussianModel()
{
    m_GaussianModel.clear();
    m_Alphas.clear();

    for ( unsigned int i = 0; i < m_NbClasses; i++ )
    {
        GaussianFunctionType::Pointer densityFunction1 = GaussianFunctionType::New();
        densityFunction1->SetMeasurementVectorSize( m_NbClasses );
        GaussianFunctionType::MeanVectorType mean1( m_NbClasses );
        GaussianFunctionType::CovarianceMatrixType cov1;
        cov1.SetSize( m_NbClasses, m_NbClasses );
        cov1.Fill(0);

        // Take information from T1-w
        GaussianFunctionType::MeanVectorType t1mean = (m_Solution1D[i])->GetMean();
        GaussianFunctionType::CovarianceMatrixType t1var = (m_Solution1D[i])->GetCovariance();

        mean1[0] = t1mean[0];
        mean1[1] = m_SelectedMax[1][i];
        mean1[2] = m_SelectedMax[2][i];

        cov1[0][0] = t1var[0][0];
        cov1[1][1] = m_Stds[1][i]*m_Stds[1][i];
        cov1[2][2] = m_Stds[2][i]*m_Stds[2][i];

        densityFunction1->SetMean( mean1 );
        densityFunction1->SetCovariance( cov1 );

        m_GaussianModel.push_back(densityFunction1);
        m_Alphas.push_back(m_Solution1DAlphas[i]);
    }
}

template <typename TInputImage, typename TMaskImage>
float HierarchicalInitializer<TInputImage,TMaskImage>
::regionMedianValue(itk::Image< float, 3 >::Pointer image, typename TInputImage::Pointer mask )
{
    ImageIteratorTypeF imageIt(image, image->GetLargestPossibleRegion() );
    InputIteratorType maskIt(mask, mask->GetLargestPossibleRegion() );

    std::vector<float> regionValueVector;

    while(!imageIt.IsAtEnd())
    {
        if ( maskIt.Get() != 0 )
        {
            regionValueVector.push_back(imageIt.Get());
        }
        ++imageIt;
        ++maskIt;
    }
    std::sort (regionValueVector.begin(), regionValueVector.end());
    return regionValueVector[regionValueVector.size()/2];
}

template <typename TInputImage, typename TMaskImage>
void HierarchicalInitializer<TInputImage,TMaskImage>::Update()
{
    m_ImagesVector.clear();

    if(m_IndexImage1<m_NbMaxImages){m_ImagesVector.push_back(this->GetInputImage1());}
    if(m_IndexImage2<m_NbMaxImages){m_ImagesVector.push_back(this->GetInputImage2());}
    if(m_IndexImage3<m_NbMaxImages){m_ImagesVector.push_back(this->GetInputImage3());}
    if(m_IndexImage4<m_NbMaxImages){m_ImagesVector.push_back(this->GetInputImage4());}
    if(m_IndexImage5<m_NbMaxImages){m_ImagesVector.push_back(this->GetInputImage5());}

    m_NumberOfModalities = m_ImagesVector.size();

    this->ComputeSolution1D();
    this->ComputeInitialT1Classification();
    this->ComputeHistogram();
    this->ComputeVariances();
    this->FillInitialGaussianModel();
}


} // end of namespace
