#pragma once

#include "animaNLinksFilter.h"

namespace anima
{

template <typename TInput, typename TOutput>
void NLinksFilter<TInput, TOutput>::SetInputImage1(const TInput* image)
{
    this->SetNthInput(0, const_cast<TInput*>(image));
    m_NbModalities++;
    m_ListImages.push_back(image);
}

template <typename TInput, typename TOutput>
void NLinksFilter<TInput, TOutput>::SetInputSeedProbaSources(const TSeedProba* proba)
{
    this->SetNthInput(1, const_cast<TSeedProba*>(proba));
}

template <typename TInput, typename TOutput>
void NLinksFilter<TInput, TOutput>::SetInputSeedProbaSinks(const TSeedProba* proba)
{
    this->SetNthInput(2, const_cast<TSeedProba*>(proba));
}

template <typename TInput, typename TOutput>
void NLinksFilter<TInput, TOutput>::SetMask(const TMask* mask)
{
    this->SetNthInput(3, const_cast<TMask*>(mask));
}

template <typename TInput, typename TOutput>
void NLinksFilter<TInput, TOutput>::SetInputImage2(const TInput* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInput*>(image));
    m_IndexImage2 = m_NbInputs;
    m_NbInputs++;
    m_NbModalities++;
    m_ListImages.push_back(image);
}

template <typename TInput, typename TOutput>
void NLinksFilter<TInput, TOutput>::SetInputImage3(const TInput* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInput*>(image));
    m_IndexImage3 = m_NbInputs;
    m_NbInputs++;
    m_NbModalities++;
    m_ListImages.push_back(image);
}

template <typename TInput, typename TOutput>
void NLinksFilter<TInput, TOutput>::SetInputImage4(const TInput* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInput*>(image));
    m_IndexImage4 = m_NbInputs;
    m_NbInputs++;
    m_NbModalities++;
    m_ListImages.push_back(image);
}

template <typename TInput, typename TOutput>
void NLinksFilter<TInput, TOutput>::SetInputImage5(const TInput* image)
{
    this->SetNthInput(m_NbInputs, const_cast<TInput*>(image));
    m_IndexImage5 = m_NbInputs;
    m_NbInputs++;
    m_NbModalities++;
    m_ListImages.push_back(image);
}


template <typename TInput, typename TOutput>
typename TInput::ConstPointer NLinksFilter<TInput, TOutput>::GetInputImage1()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(0) );
}

template <typename TInput, typename TOutput>
itk::Image <float,3>::ConstPointer NLinksFilter<TInput, TOutput>::GetInputSeedProbaSources()
{
    return static_cast< const TSeedProba * >
            ( this->itk::ProcessObject::GetInput(1) );
}

template <typename TInput, typename TOutput>
itk::Image <float,3>::ConstPointer NLinksFilter<TInput, TOutput>::GetInputSeedProbaSinks()
{
    return static_cast< const TSeedProba * >
            ( this->itk::ProcessObject::GetInput(2) );
}

template <typename TInput, typename TOutput>
itk::Image <unsigned char,3>::ConstPointer NLinksFilter<TInput, TOutput>::GetMask()
{
    return static_cast< const TMask * >
            ( this->itk::ProcessObject::GetInput(3) );
}

template <typename TInput, typename TOutput>
typename TInput::ConstPointer NLinksFilter<TInput, TOutput>::GetInputImage2()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage2) );
}

template <typename TInput, typename TOutput>
typename TInput::ConstPointer NLinksFilter<TInput, TOutput>::GetInputImage3()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage3) );
}

template <typename TInput, typename TOutput>
typename TInput::ConstPointer NLinksFilter<TInput, TOutput>::GetInputImage4()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage4) );
}

template <typename TInput, typename TOutput>
typename TInput::ConstPointer NLinksFilter<TInput, TOutput>::GetInputImage5()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(m_IndexImage5) );
}


template <typename TInputImage, typename TOutput>
itk::DataObject::Pointer NLinksFilter <TInputImage, TOutput>::MakeOutput(unsigned int idx)
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
typename TOutput::Pointer NLinksFilter <TInputImage, TOutput>::GetOutput()
{
    return dynamic_cast< TOutput* >( this->itk::ProcessObject::GetOutput(0) );
}
template <typename TInputImage, typename TOutput>
typename TOutput::Pointer NLinksFilter <TInputImage, TOutput>::GetOutputBackground()
{
    return dynamic_cast< TOutput* >( this->itk::ProcessObject::GetOutput(1) );
}

template <typename TInput, typename TOutput>
bool NLinksFilter<TInput, TOutput>::readMatrixFile()
{
    // Read and Parse the data
    typedef itk::CSVArray2DFileReader<float> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName ( m_MatFilename );
    reader->SetFieldDelimiterCharacter( ';' );
    reader->SetStringDelimiterCharacter( '"' );
    reader->HasColumnHeadersOff();
    reader->HasRowHeadersOff();

    // Try Non-existent filename
    try
    {
        reader->Parse();
    }
    catch( itk::ExceptionObject & exp )
    {
        std::cerr << "Expected Exception caught!" << std::endl;
        std::cerr << exp << std::endl;
        return false;
    }

    reader->Print(std::cout);

    typedef itk::CSVArray2DDataObject<float> DataFrameObjectType;
    DataFrameObjectType::Pointer dfo = reader->GetOutput();

    m_Matrix = dfo->GetMatrix();
    return true;
}

template <typename TInput, typename TOutput>
void
NLinksFilter<TInput, TOutput>
::CheckSpectralGradient()
{
    unsigned int nbInputs = this->GetNumberOfIndexedInputs();
    if (nbInputs < 4)
    {
        std::cerr << "Error: No inputs available... Exiting..." << std::endl;
        exit(-1);
    }

    if((m_MatFilename!="") && m_UseSpectralGradient)
    {
        if(m_ListImages.size()<3)
        {
            std::cerr << "-- Error in Graph3D: Wrong number of images, must be equal or superior to 3" <<  std::endl;
            exit(-1);
        }
        if(!readMatrixFile())
        {
            std::cerr << "-- Error in Graph3D: error while reading matrix file" << std::endl;
            exit(-1);
        }
        if( (m_Matrix.Rows()!=3) || (m_Matrix.Cols()!=m_ListImages.size()) )
        {
            std::cerr << "-- Error in Graph3D: Wrong number of coefficients in matrix after reading matrix csv file" <<  std::endl;
            exit(-1);
        }
        std::cout << "using matfilename" << std::endl;
    }
    else if( (m_Matrix.Rows()!=0) && m_UseSpectralGradient)
    {
        if( (m_Matrix.Rows()!=3) || (m_Matrix.Cols()!=m_ListImages.size()) )
        {
            std::cerr << "-- Error in NLinksFilter: Wrong number of coefficients in spectral gradient matrix" <<  std::endl;
            std::cerr << "-- Matrix.Rows(): " << m_Matrix.Rows() << std::endl;
            std::cerr << "-- Matrix.Cols(): " << m_Matrix.Cols() << std::endl;
            std::cerr << "-- ListImages.size(): " << m_ListImages.size()<<  std::endl;
            exit(-1);
        }

        std::cout << "using given matrix" << std::endl;
        std::cout << m_Matrix(0,0) << " "<< m_Matrix(0,1) << " "<< m_Matrix(0,2) << std::endl;
        std::cout << m_Matrix(1,0) << " "<< m_Matrix(1,1) << " "<< m_Matrix(1,2) << std::endl;
        std::cout << m_Matrix(2,0) << " "<< m_Matrix(2,1) << " "<< m_Matrix(2,2) << std::endl;
    }
    else
    {
        // If no matrix is given, set orginal spectral gradient
        m_Matrix.SetSize(3,3);

        m_Matrix(0,0) = 0.002358; m_Matrix(0,1) = 0.025174;  m_Matrix(0,2) = 0.010821;
        m_Matrix(1,0) = 0.011943; m_Matrix(1,1) = 0.001715;  m_Matrix(1,2) = -0.013994;
        m_Matrix(2,0) = 0.013743; m_Matrix(2,1) = -0.023965; m_Matrix(2,2) = 0.00657;

        // Check if possible use of spectral gradient
        if (m_UseSpectralGradient && (m_ListImages.size() != 3))
        {
            std::cerr << "-- Warning in Graph3D: the spectral gradient requires 3 modalities" << std::endl << "Falling back to conventional gradient approach"  << std::endl;
            m_UseSpectralGradient = false;
        }

        if (m_UseSpectralGradient)
            std::cout << "using original grad spec" << std::endl;
        else
            std::cout << "NO using grad spec" << std::endl;
    }

    m_size = this->GetMask()->GetLargestPossibleRegion().GetSize();

}

template <typename TInput, typename TOutput>
void
NLinksFilter<TInput, TOutput>
::GenerateData()
{
    typename TOutput::Pointer output = this->GetOutput();
    output->SetRegions(this->GetMask()->GetLargestPossibleRegion());
    output->CopyInformation(this->GetMask());
    output->Allocate();
    output->FillBuffer(0);

    typename TOutput::Pointer outputBackground = this->GetOutputBackground();
    outputBackground->SetRegions(this->GetMask()->GetLargestPossibleRegion());
    outputBackground->CopyInformation(this->GetMask());
    outputBackground->Allocate();
    outputBackground->FillBuffer(0);

    std::cout << "Computing N-Links..."<< std::endl;

    this->CheckSpectralGradient();
    this->CreateGraph();
    this->SetGraph();

    m_graph -> maxflow();

    int cpt=0;
    MaskRegionConstIteratorType maskIt (this->GetMask(),this->GetMask()->GetLargestPossibleRegion());
    OutputIteratorType outIt (output,output->GetLargestPossibleRegion() );
    OutputIteratorType outBackgroundIt (outputBackground,outputBackground->GetLargestPossibleRegion() );

    while (!maskIt.IsAtEnd())
    {
        outIt.Set(0);
        if (maskIt.Get() != 0)
        {
            int v = m_graph->what_segment(cpt);
            unsigned char buff = 0 | (v == GraphType::SOURCE) ? 1 : 0;
            outIt.Set(static_cast<OutputPixelType>(buff));
            outBackgroundIt.Set(1-buff);
            cpt++;
        }
        ++maskIt;
        ++outBackgroundIt;
        ++outIt;
    }

    m_NbModalities = 0;
    m_NbInputs = 3;
    m_ListImages.clear();
    if (m_graph) delete m_graph;

}


template <typename TInput, typename TOutput>
double NLinksFilter<TInput, TOutput>::computeNLink(int i1, int j1, int k1, int i2, int j2, int k2)
{
    pixelIndexInt index1, index2;
    index1[0]=i1;
    index1[1]=j1;
    index1[2]=k1;

    index2[0]=i2;
    index2[1]=j2;
    index2[2]=k2;

    if (m_UseSpectralGradient)
    {
        double g = m_e1->GetPixel(index1) - m_e1->GetPixel(index2);
        double w = m_e2->GetPixel(index1) - m_e2->GetPixel(index2);
        return static_cast<double>((.1+std::exp(-(g * g + w * w) / (2 * m_Sigma * m_Sigma))));
    }
    else
    {
        int dim = m_ListImages.size();
        double g;
        double sum = 0.0;
        for (int m=0; m<dim; m++)
        {
            g = m_ListImages[m]->GetPixel(index1) - m_ListImages[m]->GetPixel(index2);
            sum += g * g;
        }
        return static_cast<double>((.1+std::exp(- sum / (2 * m_Sigma * m_Sigma))));
    }
}

template <typename TInput, typename TOutput>
void NLinksFilter<TInput, TOutput>::CreateGraph()
{
    // Compute e, e_l, e_l_l
    // These quantities are useful to compute the spectral gradient
    if(m_UseSpectralGradient)
    {
        m_e1 = TSeedProba::New();
        m_e1->SetRegions(this->GetMask()->GetLargestPossibleRegion());
        m_e1->CopyInformation(this->GetMask());
        m_e1->Allocate();
        m_e1->FillBuffer(0);

        m_e2 = TSeedProba::New();
        m_e2->SetRegions(this->GetMask()->GetLargestPossibleRegion());
        m_e2->CopyInformation(this->GetMask());
        m_e2->Allocate();
        m_e2->FillBuffer(0);

        ImageIteratorTypeProba e1It (m_e1, m_e1->GetLargestPossibleRegion() );
        ImageIteratorTypeProba e2It (m_e2, m_e2->GetLargestPossibleRegion() );

        std::vector<InRegionConstIteratorType > ListImagesVectorIt;
        for ( unsigned int i = 0; i < m_ListImages.size(); i++ )
        {
            InRegionConstIteratorType It(m_ListImages[i],m_ListImages[i]->GetLargestPossibleRegion() );
            ListImagesVectorIt.push_back(It);
        }

        while((!e1It.IsAtEnd()) && (!e2It.IsAtEnd()))
        {

            float e = 0;
            float e_l = 0;
            float e_l_l = 0;

            for(unsigned int m=0; m < m_ListImages.size(); m++)
            {
                e +=     m_Matrix(0,m) * (ListImagesVectorIt[m].Get());
                e_l +=   m_Matrix(1,m) * (ListImagesVectorIt[m].Get());
                e_l_l += m_Matrix(2,m) * (ListImagesVectorIt[m].Get());
            }
            e1It.Set(e_l/e);
            e2It.Set((e*e_l_l-std::pow(e_l, 2.f)) / std::pow(e,2.f) );

            ++e1It;
            ++e2It;
            for ( unsigned int i = 0; i < m_ListImages.size(); i++ )
            {
                ++ListImagesVectorIt[i];
            }
        }

    }
}

template <typename TInput, typename TOutput>
bool NLinksFilter<TInput, TOutput>::isInside (unsigned int x, unsigned int y, unsigned int z ) const
{
    return ( x < m_size[0] && y < m_size[1] && z < m_size[2] &&  x >= 0 && y >= 0 && z >= 0 );
}

template <typename TInput, typename TOutput>
void NLinksFilter<TInput, TOutput>::SetGraph()
{
    // allocate only necessary memory
    int nb_vox = 0;

    MaskRegionConstIteratorType maskIt (this->GetMask(),this->GetMask()->GetLargestPossibleRegion() );
    while (!maskIt.IsAtEnd())
    {
        if (maskIt.Get() != 0)
        {
            nb_vox++;
        }
        ++maskIt;
    }

    int nb_edges = 7*nb_vox;

    try
    {
        m_graph = new GraphType(nb_vox, nb_edges);
    }
    catch (std::bad_alloc& ba)
    {
        std::cerr << "-- Error in NLinksFilter: insufficient memory to create the graph: " << ba.what() << '\n';
        exit(-1);
    }

    // Create the nodes of the graph and set the correspondences with the original image
    int compt = 0;
    ImageTypeInt::Pointer pix = ImageTypeInt::New();
    pix->SetRegions(this->GetMask()->GetLargestPossibleRegion());
    pix->CopyInformation(this->GetMask());
    pix->Allocate();
    pix->FillBuffer(-1);

    ImageIteratorTypeInt pixIt (pix,pix->GetLargestPossibleRegion());

    maskIt.GoToBegin();
    while (!maskIt.IsAtEnd())
    {
        if (maskIt.Get() != 0)
        {
            m_graph -> add_node();
            pixIt.Set(compt++);
        }
        ++maskIt;
        ++pixIt;
    }


    maskIt.GoToBegin();
    pixIt.GoToBegin();


    // Create the t-links and n-links
    while (!maskIt.IsAtEnd())
    {
        if(maskIt.Get() != 0)
        {
            pixelIndexInt index, index1;
            index[0]=maskIt.GetIndex()[0];
            index[1]=maskIt.GetIndex()[1];
            index[2]=maskIt.GetIndex()[2];

            int pix_ref = pix->GetPixel(index);
            int pix_courant;
            double cap, rcap;

            // Compute the 6 n-links of each standard node (gradients between the current voxel and its neighbors)
            index1[0]=index[0]+1;
            index1[1]=index[1];
            index1[2]=index[2];
            if ( (isInside(index[0]+1, index[1] , index[2])) && (this->GetMask()->GetPixel(index1)!=0) )
            {
                pix_courant = pix->GetPixel(index1);
                cap  = computeNLink(index[0]+1, index[1], index[2], index[0], index[1], index[2]);
                rcap = cap;
                if (!(cap >= 0))
                    cap = rcap = 0;
                m_graph -> add_edge(pix_ref, pix_courant, cap, rcap);
            }


            index1[0]=index[0];
            index1[1]=index[1]+1;
            index1[2]=index[2];
            if ( (isInside(index[0], index[1]+1 , index[2])) && (this->GetMask()->GetPixel(index1)!=0) )
            {
                pix_courant = pix->GetPixel(index1);
                cap  = computeNLink(index[0], index[1]+1, index[2], index[0], index[1], index[2]);
                rcap = cap;
                if (!(cap >= 0)) // eq cap < 0 ? no ???
                    cap = rcap = 0;
                m_graph -> add_edge(pix_ref, pix_courant, cap, rcap);
            }

            index1[0]=index[0];
            index1[1]=index[1];
            index1[2]=index[2]+1;
            if ( (isInside(index[0], index[1] , index[2]+1)) && (this->GetMask()->GetPixel(index1)!=0) )
            {
                pix_courant = pix->GetPixel(index1);
                cap  = computeNLink(index[0], index[1], index[2]+1, index[0], index[1], index[2]);
                rcap = cap;
                if (!(cap >= 0))
                    cap = rcap = 0;
                m_graph -> add_edge(pix_ref, pix_courant, cap, rcap);
            }

            // Create the t-links to the source and the sink
            double t_source = static_cast<double>(this->GetInputSeedProbaSources()->GetPixel(index));
            double t_sink   = static_cast<double>(this->GetInputSeedProbaSinks()->GetPixel(index));
            m_graph -> add_tweights(pix_ref, t_source, t_sink);

        }

        ++pixIt;
        ++maskIt;
    }
}

} //end of namespace anima
