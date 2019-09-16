#pragma once

#include "animaComputeSolution.h"

namespace anima
{

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
void ComputeSolution <TInputImage, TMaskImage, TAtlasImage>::SetMask(const TMaskImage* mask)
{
    this->SetNthInput(0, const_cast<TMaskImage*>(mask));
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
void ComputeSolution <TInputImage, TMaskImage, TAtlasImage>::SetInputImage1(const TInputImage* image)
{
    this->SetNthInput(1, const_cast<TInputImage*>(image));
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
void ComputeSolution <TInputImage, TMaskImage, TAtlasImage>::SetInputImage2(const TInputImage* image)
{
    this->SetNthInput(2, const_cast<TInputImage*>(image));
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
void ComputeSolution <TInputImage, TMaskImage, TAtlasImage>::SetInputImage3(const TInputImage* image)
{
    this->SetNthInput(3, const_cast<TInputImage*>(image));
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
void ComputeSolution <TInputImage, TMaskImage, TAtlasImage>::SetInputCSFAtlas(const TAtlasImage* atlas)
{
    this->SetNthInput(4, const_cast<TAtlasImage*>(atlas));
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
void ComputeSolution <TInputImage, TMaskImage, TAtlasImage>::SetInputGMAtlas(const TAtlasImage* atlas)
{
    this->SetNthInput(5, const_cast<TAtlasImage*>(atlas));
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
void ComputeSolution <TInputImage, TMaskImage, TAtlasImage>::SetInputWMAtlas(const TAtlasImage* atlas)
{
    this->SetNthInput(6, const_cast<TAtlasImage*>(atlas));
}



template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
typename TMaskImage::ConstPointer ComputeSolution <TInputImage, TMaskImage, TAtlasImage>::GetMask()
{
    return static_cast< const ImageTypeUC * >
            ( this->itk::ProcessObject::GetInput(0) );
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
typename TInputImage::ConstPointer ComputeSolution <TInputImage, TMaskImage, TAtlasImage>::GetInputImage1()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(1) );
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
typename TInputImage::ConstPointer ComputeSolution <TInputImage, TMaskImage, TAtlasImage>::GetInputImage2()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(2) );
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
typename TInputImage::ConstPointer ComputeSolution <TInputImage, TMaskImage, TAtlasImage>::GetInputImage3()
{
    return static_cast< const TInputImage * >
            ( this->itk::ProcessObject::GetInput(3) );
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
typename TAtlasImage::ConstPointer ComputeSolution <TInputImage, TMaskImage, TAtlasImage>::GetInputCSFAtlas()
{
    return static_cast< const TAtlasImage * >
            ( this->itk::ProcessObject::GetInput(4) );
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
typename TAtlasImage::ConstPointer ComputeSolution <TInputImage, TMaskImage, TAtlasImage>::GetInputGMAtlas()
{
    return static_cast< const TAtlasImage * >
            ( this->itk::ProcessObject::GetInput(5) );
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
typename TAtlasImage::ConstPointer ComputeSolution <TInputImage, TMaskImage, TAtlasImage>::GetInputWMAtlas()
{
    return static_cast< const TAtlasImage * >
            ( this->itk::ProcessObject::GetInput(6) );
}


template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
void
ComputeSolution <TInputImage, TMaskImage, TAtlasImage>::SortGaussianModel()
{
    std::vector<GaussianFunctionType::Pointer> tmpGaussianModel;
    std::vector<double> tmpAlphas;

    typedef std::map<double,unsigned int> OrderType;
    OrderType ordered;
    OrderType::iterator it;

    ordered.clear();
    for(unsigned int i = 0; i < m_GaussianModel.size(); i++)
    {
        GaussianFunctionType::MeanVectorType mean = (m_GaussianModel[i])->GetMean();
        ordered.insert(OrderType::value_type(mean[0],i));
    }

    for(it = ordered.begin(); it != ordered.end(); ++it)
    {
        tmpGaussianModel.push_back(m_GaussianModel[it->second]);
        tmpAlphas.push_back(m_Alphas[it->second]);
    }
    m_GaussianModel = tmpGaussianModel;
    m_Alphas = tmpAlphas;
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
void
ComputeSolution <TInputImage, TMaskImage, TAtlasImage>
::CheckInputs()
{
    if( (this->GetMask().IsNull()) || (this->GetInputImage1().IsNull()) || (this->GetInputImage2().IsNull()) || (this->GetInputImage3().IsNull()) )
    {
        std::cerr << "Error: Inputs are missing... Exiting..." << std::endl;
        exit(-1);
    }

    if(m_UseT2 && m_UseDP && m_UseFLAIR)
    {
        std::cerr << "-- Error in Automatic segmentation: only 2 images among T2, DP and FLAIR must be used for automatic segmentation" << std::endl;
        exit(-1);
    }

    if(m_InitMethodType==0)
    {
        if((this->GetInputCSFAtlas().IsNull()) || (this->GetInputGMAtlas().IsNull()) || (this->GetInputWMAtlas().IsNull()))
        {
            std::cerr << "Error: Some atlas images are missing for the initialization... Exiting..." << std::endl;
            exit(-1);
        }
    }

    if(m_InitMethodType==1)
    {
        if((m_UseT2==false) || (m_UseDP==false))
        {
            std::cerr << "Error: Automatic segmentation with Hierarchical DP initialisation requires T2 and DP images... Exiting..." << std::endl;
            exit(-1);
        }
    }

    if(m_InitMethodType==2)
    {
        if(m_UseFLAIR==false)
        {
            std::cerr << "Error: Automatic segmentation with Hierarchical FLAIR initialisation requires FLAIR image... Exiting..." << std::endl;
            exit(-1);
        }
    }
}


template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
void
ComputeSolution <TInputImage, TMaskImage, TAtlasImage>
::RescaleImages()
{
    m_InputImage_T1_UC = ImageTypeUC::New();
    m_InputImage_T2_DP_UC = ImageTypeUC::New();
    m_InputImage_DP_FLAIR_UC = ImageTypeUC::New();

    double desiredMinimum=0,desiredMaximum=255;

    typename RescaleFilterType::Pointer rescaleFilter1 = RescaleFilterType::New();
    rescaleFilter1->SetInput( this->GetInputImage1() );
    rescaleFilter1->SetOutputMinimum( desiredMinimum );
    rescaleFilter1->SetOutputMaximum( desiredMaximum );
    rescaleFilter1->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    rescaleFilter1->UpdateLargestPossibleRegion();
    m_InputImage_T1_UC = rescaleFilter1->GetOutput();

    typename RescaleFilterType::Pointer rescaleFilter2 = RescaleFilterType::New();
    rescaleFilter2->SetInput( this->GetInputImage2() );
    rescaleFilter2->SetOutputMinimum( desiredMinimum );
    rescaleFilter2->SetOutputMaximum( desiredMaximum );
    rescaleFilter2->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    rescaleFilter2->UpdateLargestPossibleRegion();
    m_InputImage_T2_DP_UC = rescaleFilter2->GetOutput();

    typename RescaleFilterType::Pointer rescaleFilter3 = RescaleFilterType::New();
    rescaleFilter3->SetInput( this->GetInputImage3() );
    rescaleFilter3->SetOutputMinimum( desiredMinimum );
    rescaleFilter3->SetOutputMaximum( desiredMaximum );
    rescaleFilter3->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    rescaleFilter3->UpdateLargestPossibleRegion();
    m_InputImage_DP_FLAIR_UC = rescaleFilter3->GetOutput();
}


template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
int
ComputeSolution <TInputImage, TMaskImage, TAtlasImage>
::WriteSolution(std::string filename)
{
    const unsigned int ARows = 1;
    const unsigned int ACols = 39; // m_NbTissus*(1+m_NbModalities+m_NbModalities*m_NbModalities);

    typedef itk::Array2D<double> MatrixType;
    MatrixType matrix(ARows,ACols);

    unsigned int t = 0;
    for(unsigned int i = 0; i < m_NbTissus; i++)
    {
        matrix[0][t] = m_Alphas[i];
        t++;

        GaussianFunctionType::MeanVectorType mu = (m_GaussianModel[i])->GetMean();
        for(unsigned int j = 0; j < m_NbModalities; j++)
        {
            matrix[0][t] = mu[j];
            t++;
        }

        GaussianFunctionType::CovarianceMatrixType covar = (m_GaussianModel[i])->GetCovariance();
        for(unsigned int k = 0; k < m_NbModalities; k++)
        {
            for(unsigned int j = 0; j < m_NbModalities; j++)
            {
                matrix[0][t] = covar(k,j);
                t++;
            }
        }
    }

    // write out the array2D object
    typedef itk::CSVNumericObjectFileWriter<double, ARows, ACols> WriterType;
    WriterType::Pointer writer = WriterType::New();

    writer->SetFieldDelimiterCharacter(';');
    writer->SetFileName( filename );
    writer->SetInput( &matrix );
    try
    {
        writer->Write();
    }
    catch (itk::ExceptionObject& exp)
    {
        std::cerr << "Exception caught!" << std::endl;
        std::cerr << exp << std::endl;
        return -1;
    }

    return 0;
}


template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
int
ComputeSolution <TInputImage, TMaskImage, TAtlasImage>
::ReadSolution(std::string filename)
{
    typedef itk::Array2D<double> MatrixType;
    typedef itk::CSVArray2DFileReader<double > ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( filename );
    reader->SetFieldDelimiterCharacter( ';' );
    reader->SetStringDelimiterCharacter( '"' );
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();

    // read the file
    try
    {
        reader->Update();
    }
    catch (itk::ExceptionObject& exp)
    {
        std::cerr << "Exception caught!" << std::endl;
        std::cerr << exp << std::endl;
        return -1;
    }

    typedef itk::CSVArray2DDataObject<double> DataFrameObjectType;
    DataFrameObjectType::Pointer dfo = reader->GetOutput();
    MatrixType matrix = dfo->GetMatrix();

    unsigned int nbCols = m_NbTissus*(1+m_NbModalities+m_NbModalities*m_NbModalities);
    if(matrix.rows()!=1 || matrix.cols()!=nbCols)
    {
        std::cout<< "wrong type of matrix file... cannot read solution" << std::endl;
        return -1;
    }

    unsigned int t = 0;
    for(unsigned int i = 0; i < m_NbTissus; i++)
    {
        m_Alphas.push_back(matrix[0][t]);
        t++;

        GaussianFunctionType::Pointer densityFunction = GaussianFunctionType::New();
        densityFunction->SetMeasurementVectorSize( m_NbTissus );
        GaussianFunctionType::MeanVectorType mean( m_NbModalities );
        GaussianFunctionType::CovarianceMatrixType cov;
        cov.SetSize( m_NbModalities, m_NbModalities );
        cov.Fill(0);

        for(unsigned int j = 0; j < m_NbModalities; j++)
        {
            mean[j] = matrix[0][t];
            t++;
        }

        for(unsigned int k = 0; k < m_NbModalities; k++)
        {
            for(unsigned int j = 0; j < m_NbModalities; j++)
            {
                cov[k][j] = matrix[0][t];
                t++;
            }
        }

        densityFunction->SetMean( mean );
        densityFunction->SetCovariance( cov );
        m_GaussianModel.push_back( densityFunction );
    }

    return 0;
}

template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
int
ComputeSolution <TInputImage, TMaskImage, TAtlasImage>
::PrintSolution(std::vector<double> alphas, std::vector<GaussianFunctionType::Pointer> model)
{
    unsigned int nbTissus = alphas.size();
    unsigned int nbModalities = (model[0])->GetMean().Size();
    for(unsigned int i = 0; i < nbTissus; i++)
    {
        std::cout << "* Class: " << i << std::endl;
        std::cout << "  Alpha: " << alphas[i] << std::endl;
        GaussianFunctionType::MeanVectorType mu = (model[i])->GetMean();
        std::cout << "  Mean: " << std::endl;
        for(unsigned int j = 0; j < nbModalities; j++)
        {
            std::cout << "    " << mu[j] << std::endl;
        }

        GaussianFunctionType::CovarianceMatrixType covar = (model[i])->GetCovariance();
        std::cout << "  covar: " << std::endl;
        for(unsigned int k = 0; k < nbModalities; k++)
        {
            std::cout << "    " << covar(k,0) << " " << covar(k,1) << " " << covar(k,2) << std::endl;
        }
    }
    return 0;
}


template <typename TInputImage, typename TMaskImage, typename TAtlasImage>
void
ComputeSolution <TInputImage, TMaskImage, TAtlasImage>
::Update()
{
    this->CheckInputs();
    this->RescaleImages();

    if((m_GaussianModel.size() == m_NbTissus) && (m_Alphas.size() == m_NbTissus))
    {
        m_SolutionSet = true;
        std::cout << "Solution already set..."<< std::endl;
        if( m_Verbose )
        {
            this->PrintSolution(m_Alphas, m_GaussianModel);
        }
    }
    else
    {
        if(m_SolutionReadFilename!="")
        {
            m_SolutionSet = true;
            std::cout << "Reading solution file..." << std::endl;
            if(ReadSolution(m_SolutionReadFilename))
            {
                m_SolutionSet = false;
            }
            if( m_Verbose )
            {
                this->PrintSolution(m_Alphas, m_GaussianModel);
            }
        }
    }

    if(!m_SolutionSet)
    {
        // Choose intializer
        ModelInitializer::Pointer initializer;
        switch (m_InitMethodType)
        {
        case 0:
        {
            initializer = AtlasInitializerType::New();
            dynamic_cast<AtlasInitializerType *>( initializer.GetPointer() ) ->SetMask( this->GetMask() );
            dynamic_cast<AtlasInitializerType *>( initializer.GetPointer() ) ->SetInputImage1( m_InputImage_T1_UC );
            dynamic_cast<AtlasInitializerType *>( initializer.GetPointer() ) ->SetInputImage2( m_InputImage_T2_DP_UC );
            dynamic_cast<AtlasInitializerType *>( initializer.GetPointer() ) ->SetInputImage3( m_InputImage_DP_FLAIR_UC );
            dynamic_cast<AtlasInitializerType *>( initializer.GetPointer() ) ->SetAtlasImage1( this->GetInputCSFAtlas() );
            dynamic_cast<AtlasInitializerType *>( initializer.GetPointer() ) ->SetAtlasImage2( this->GetInputGMAtlas() );
            dynamic_cast<AtlasInitializerType *>( initializer.GetPointer() ) ->SetAtlasImage3( this->GetInputWMAtlas() );
            std::cout<< "Choosen initializer: Atlas" << std::endl;
            break;
        }
        case 1:
        {
            bool use_HierarFLAIR = false;
            initializer = HierarchicalType::New();
            dynamic_cast<HierarchicalType *>( initializer.GetPointer() ) ->SetMask( this->GetMask() );
            dynamic_cast<HierarchicalType *>( initializer.GetPointer() ) ->SetInputImage1( m_InputImage_T1_UC );
            dynamic_cast<HierarchicalType *>( initializer.GetPointer() ) ->SetInputImage2( m_InputImage_T2_DP_UC );
            dynamic_cast<HierarchicalType *>( initializer.GetPointer() ) ->SetInputImage3( m_InputImage_DP_FLAIR_UC );
            dynamic_cast<HierarchicalType *>( initializer.GetPointer() ) ->SetThirdIsFLAIR( use_HierarFLAIR );
            dynamic_cast<HierarchicalType *>( initializer.GetPointer() ) ->SetRobust( m_RejRatioHierar );
            dynamic_cast<HierarchicalType *>( initializer.GetPointer() ) ->SetTol( m_Tol );
            std::cout<< "Choosen initializer: Hierarchical DP " << std::endl;
            break;
        }
        case 2:
        {
            bool use_HierarFLAIR = true;
            initializer = HierarchicalType::New();
            dynamic_cast<HierarchicalType *>( initializer.GetPointer() ) ->SetMask( this->GetMask() );
            dynamic_cast<HierarchicalType *>( initializer.GetPointer() ) ->SetInputImage1( m_InputImage_T1_UC );
            dynamic_cast<HierarchicalType *>( initializer.GetPointer() ) ->SetInputImage2( m_InputImage_T2_DP_UC );
            dynamic_cast<HierarchicalType *>( initializer.GetPointer() ) ->SetInputImage3( m_InputImage_DP_FLAIR_UC );
            dynamic_cast<HierarchicalType *>( initializer.GetPointer() ) ->SetThirdIsFLAIR( use_HierarFLAIR );
            dynamic_cast<HierarchicalType *>( initializer.GetPointer() ) ->SetRobust( m_RejRatioHierar );
            dynamic_cast<HierarchicalType *>( initializer.GetPointer() ) ->SetTol( m_Tol );
            std::cout<< "Choosen initializer: Hierarchical FLAIR" << std::endl;
            break;
        }
        default:
        {
            std::cerr<< "-- Error in compute solution filter: initialisation failed" << std::endl;
            exit(-1);
        }
        }//switch InitMethod

        std::cout << "Computing initialization for EM..." << std::endl;
        initializer->Update();
        std::vector<GaussianFunctionType::Pointer> initia = initializer->GetInitialization();
        std::vector<double> initiaAlphas = initializer->GetAlphas();

        std::cout << "Computing gaussian model..." << std::endl;
        typename GaussianREMEstimatorType::Pointer estimator = GaussianREMEstimatorType::New();
        estimator ->SetMaxIterationsConc( m_EmIter_concentration );
        estimator ->SetStremMode( m_EM_before_concentration );
        estimator ->SetRejectionRatio( m_RejRatio );
        estimator ->SetMaxIterations( m_EmIter );
        estimator ->SetModelMinDistance( m_MinDistance );
        estimator ->SetMask( this->GetMask() );
        estimator ->SetInputImage1( m_InputImage_T1_UC );
        estimator ->SetInputImage2( m_InputImage_T2_DP_UC );
        estimator ->SetInputImage3( m_InputImage_DP_FLAIR_UC );
        estimator ->SetVerbose( m_Verbose );

        itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
        callback ->SetCallback(eventCallback);
        estimator ->AddObserver(itk::ProgressEvent(), callback );
        estimator ->SetInitialGaussianModel(initia);
        estimator ->SetInitialAlphas(initiaAlphas);
        estimator ->Update();

        m_GaussianModel = estimator->GetGaussianModel();
        m_Alphas = estimator->GetAlphas();

        this->SortGaussianModel();

        if( m_Verbose )
        {
            std::cout << std::endl;
            std::cout<< "EM summary: " << std::endl << std::endl;
            std::cout<< "* Initial model: " << std::endl;
            PrintSolution(initiaAlphas, initia);
            std::cout << std::endl;
            std::cout<< "*  Final model: " << std::endl;
            this->PrintSolution(m_Alphas, m_GaussianModel);
            std::cout << std::endl;
        }
    }

    // Print NABT solution if necessary
    if(m_SolutionWriteFilename!="")
    {
        std::cout << "Writing solution in csv file..." << std::endl;
        WriteSolution(m_SolutionWriteFilename);
    }
}

} //end of namespace anima
