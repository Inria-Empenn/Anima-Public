#pragma once

#include "animaConnectedComponentsVolumeFilter.h"

namespace anima
{
template <typename ImageType>
void
ConnectedComponentVolumeFilter<ImageType>
::WriteOutputs()
{
    if( m_Outputfilename != "" )
    {
        std::cout << "Writing output image to: " << m_Outputfilename << std::endl;
        anima::writeImage<ImageType>(m_Outputfilename, this->GetOutput());
    }
}

template <typename ImageType>
void
ConnectedComponentVolumeFilter<ImageType>
::ComputeSpacing()
{
    // Compute Image Spacing, 4th dimension is not physical but temporal
    m_Spacing = this->GetInput()->GetSpacing();
    m_SpacingTot = m_Spacing[0];
    for (unsigned int i = 1; i < std::min(m_Dimension,(unsigned int)3);++i)
        m_SpacingTot *= m_Spacing[i];
}

template <typename ImageType>
void
ConnectedComponentVolumeFilter<ImageType>
::ComputeMinimumLesionSize()
{
    // Compute minimum lesion size in number of voxels
    double epsilon = 10e-6;
    double minSizeInVoxelD = m_MinSizeMM3 / m_SpacingTot;
    minSizeInVoxelD -= epsilon;
    double minSizeInVoxelD_ceil = std::ceil( minSizeInVoxelD );
    m_MinSizeVoxel = static_cast<unsigned int>( minSizeInVoxelD_ceil );
}

template <typename ImageType>
void
ConnectedComponentVolumeFilter<ImageType>
::GenerateData()
{
    this->ComputeSpacing();
    this->ComputeMinimumLesionSize();

    typename ConnectedComponentFilterType::Pointer connectedComponentFilter = ConnectedComponentFilterType::New();
    connectedComponentFilter->SetInput( this->GetInput() );
    connectedComponentFilter->SetFullyConnected( m_FullyConnected );
    if(this->GetNumberOfThreads() > 0)
        connectedComponentFilter->SetNumberOfThreads( this->GetNumberOfThreads() );
    connectedComponentFilter->SetCoordinateTolerance(m_Tol);
    connectedComponentFilter->SetDirectionTolerance(m_Tol);

    typename RelabelComponentFilterType::Pointer relabelFilter = RelabelComponentFilterType::New();
    relabelFilter->SetInput( connectedComponentFilter->GetOutput() );
    relabelFilter->SetMinimumObjectSize( m_MinSizeVoxel );
    if(this->GetNumberOfThreads() > 0)
        relabelFilter->SetNumberOfThreads( this->GetNumberOfThreads() );
    relabelFilter->SetCoordinateTolerance(m_Tol);
    relabelFilter->SetDirectionTolerance(m_Tol);
    relabelFilter->GraftOutput( this->GetOutput() );
    relabelFilter->Update();
    this->GraftOutput( relabelFilter->GetOutput() );

    m_OriginalNumberOfObjects = relabelFilter->GetOriginalNumberOfObjects();
    m_NumberOfObjects = relabelFilter->GetNumberOfObjects();

    std::vector<unsigned int> voxel_volume_lesion(relabelFilter->GetNumberOfObjects()+1, 0);
    m_vect_volume_lesionMM3.resize(relabelFilter->GetNumberOfObjects()+1, 0);
    m_TotalVolume = 0;

    ImageIteratorType relabelImageIt(relabelFilter->GetOutput(),relabelFilter->GetOutput()->GetLargestPossibleRegion());
    while(!relabelImageIt.IsAtEnd())
    {
        if(relabelImageIt.Get()!=0)
        {
            voxel_volume_lesion[relabelImageIt.Get()]++;
        }
        ++relabelImageIt;
    }

    for(unsigned int i = 1; i < voxel_volume_lesion.size(); i++)
    {
        m_vect_volume_lesionMM3[i] = voxel_volume_lesion[i] * m_SpacingTot;
        m_TotalVolume += m_vect_volume_lesionMM3[i];
    }

    if(m_OutputVolumefilename.size()!=0)
    {
        this->PrintVolumeFile();
    }

    this->PrintInfo();
}


template <typename ImageType>
void
ConnectedComponentVolumeFilter<ImageType>
::PrintVolumeFile()
{
    if(m_OutputVolumefilename.size()!=0)
    {
        std::ofstream fp(m_OutputVolumefilename.c_str(), std::ios_base::out | std::ios_base::trunc);
        if(!fp)
        {
            std::cerr << "cannot open output file" << m_OutputVolumefilename << std::endl;
        }
        else
        {
            fp << "Original number of objects: " << m_OriginalNumberOfObjects <<std::endl;
            fp << "Number of objects after cleaning small ones: " << m_NumberOfObjects << std::endl;
            fp << "Minimum connected component size: " << m_MinSizeMM3 << " mm3" << std::endl;
            if(m_Dimension==2)
            {
                fp << "Image Spacing: " << m_Spacing[0] << "*" << m_Spacing[1] << std::endl;
            }
            else
            {
                fp << "Image Spacing: " << m_Spacing[0] << "*" << m_Spacing[1] << "*"<< m_Spacing[2] << std::endl;
            }
            fp << "Total image spacing: " << m_SpacingTot << std::endl;

            fp << "Total volume: " << m_TotalVolume << " mm3" << std::endl;
            fp << "Volume of each connected component in (mm3): " << std::endl;
            for(unsigned int k = 1; k < m_vect_volume_lesionMM3.size(); k++)
            {
                fp << "-- Volume " << k << ": " << m_vect_volume_lesionMM3[k] << std::endl;
            }
            fp << std::endl;
            fp << std::endl;
            fp.close();
        }
    }
}

template <typename ImageType>
void
ConnectedComponentVolumeFilter<ImageType>
::PrintInfo()
{
    if( m_Verbose )
    {
        std::cout << "Original number of objects: " << m_OriginalNumberOfObjects <<std::endl;
        std::cout << "Number of objects after cleaning small ones: " << m_NumberOfObjects << std::endl;
        std::cout << "Minimum connected component size: " << m_MinSizeMM3 << " mm3" << std::endl;
        if(m_Dimension==2)
        {
            std::cout << "Image Spacing: " << m_Spacing[0] << "*" << m_Spacing[1] << std::endl;
        }
        else
        {
            std::cout << "Image Spacing: " << m_Spacing[0] << "*" << m_Spacing[1] << "*"<< m_Spacing[2] << std::endl;
        }
        std::cout << "Total image spacing: " << m_SpacingTot << std::endl;

        std::cout << "Total volume: " << m_TotalVolume << " mm3" << std::endl;

        std::cout << "Volume of each connected component in (mm3): " << std::endl;
        for(unsigned int k = 1; k < m_vect_volume_lesionMM3.size(); k++)
        {
            std::cout << "-- Volume " << k << ": " << m_vect_volume_lesionMM3[k] << std::endl;
        }
        std::cout << std::endl;
    }
}

} //end of namespace anima
