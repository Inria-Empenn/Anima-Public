#pragma once

#include "animaCheckStructureNeighborFilter.h"

namespace anima
{

template<typename TInput, typename TMask, typename TOutput>
void CheckStructureNeighborFilter< TInput, TMask, TOutput >::SetInputClassification(const TInput* image)
{
    this->SetNthInput(0, const_cast<TInput *>(image));
}

template<typename TInput, typename TMask, typename TOutput>
void CheckStructureNeighborFilter< TInput, TMask, TOutput >::SetInputMap(const TMask* image)
{
    this->SetNthInput(1, const_cast<TMask *>(image));
}

template<typename TInput, typename TMask, typename TOutput>
typename TInput::ConstPointer CheckStructureNeighborFilter< TInput, TMask, TOutput >::GetInputClassification()
{
  return static_cast< const TInput * >
         ( this->itk::ProcessObject::GetInput(0) );
}
template<typename TInput, typename TMask, typename TOutput>
typename TMask::ConstPointer CheckStructureNeighborFilter< TInput, TMask, TOutput >::GetInputMap()
{
  return static_cast< const  TMask * >
         ( this->itk::ProcessObject::GetInput(1) );
}

template<typename TInput, typename TMask, typename TOutput>
void
CheckStructureNeighborFilter< TInput, TMask, TOutput >
::WriteOutputs()
{
    if( m_OutputFilename != "" )
    {
        std::cout << "Writing lesions output image to: " << m_OutputFilename << std::endl;
        anima::writeImage<TOutput>(m_OutputFilename, this->GetOutput());
    }
}

template<typename TInput, typename TMask, typename TOutput >
void
CheckStructureNeighborFilter< TInput, TMask, TOutput >
::GenerateData()
{
    bool fullyConnected = false;

    // Create labeled map
    typename ConnectedComponentFilterType::Pointer ConnectedComponentFilter = ConnectedComponentFilterType::New();
    ConnectedComponentFilter->SetInput( this->GetInputClassification() );
    ConnectedComponentFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    ConnectedComponentFilter->SetCoordinateTolerance(m_Tol);
    ConnectedComponentFilter->SetDirectionTolerance(m_Tol);
    ConnectedComponentFilter->SetFullyConnected( fullyConnected );

    unsigned int radius = 1;

    StructuringElementType structuringElement;
    structuringElement.SetRadius(radius);
    structuringElement.CreateStructuringElement();

    DilateFilterType::Pointer dilateFilter = DilateFilterType::New();
    dilateFilter->SetInput(ConnectedComponentFilter->GetOutput());
    dilateFilter->SetKernel(structuringElement);
    dilateFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    dilateFilter->SetCoordinateTolerance(m_Tol);
    dilateFilter->SetDirectionTolerance(m_Tol);
    dilateFilter->Update();

    // Find contours
    typename LabelContourFilterType::Pointer filterLabelContour = LabelContourFilterType::New();
    filterLabelContour->SetInput(dilateFilter->GetOutput());
    filterLabelContour->SetFullyConnected( fullyConnected );
    filterLabelContour->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    filterLabelContour->SetCoordinateTolerance(m_Tol);
    filterLabelContour->SetDirectionTolerance(m_Tol);
    filterLabelContour->Update();

    unsigned int objectCount = ConnectedComponentFilter->GetObjectCount();

    std::vector<unsigned int> nbVoxelContour(objectCount+1,0);
    std::vector<unsigned int> nbVoxelIntersec(objectCount+1,0);
    std::vector<bool> enoughContour(objectCount+1,true);

    ImageIteratorTypeInt labelContourIt(filterLabelContour->GetOutput(), filterLabelContour->GetOutput()->GetLargestPossibleRegion());
    MaskConstIteratorType mapIt(this->GetInputMap(), this->GetInputMap()->GetLargestPossibleRegion());
    while(!labelContourIt.IsAtEnd())
    {
        if(labelContourIt.Get()!=0)
        {
            nbVoxelContour[labelContourIt.Get()]++;
            if(mapIt.Get()!=0)
            {
                nbVoxelIntersec[labelContourIt.Get()]++;
            }
        }
        ++labelContourIt;
        ++mapIt;
    }

    // calculer les ratios
    for(unsigned int i = 1; i < objectCount+1; i++)
    {
        if(static_cast<float>(nbVoxelIntersec[i])/static_cast<float>(nbVoxelContour[i]) < m_Ratio)
        {
            enoughContour[i] = false;
        }
    }


    // Calculer output
    OutputImagePointer output = this->GetOutput();
    output->SetRegions( this->GetInputClassification()->GetLargestPossibleRegion() );
    output->CopyInformation( this->GetInputClassification() );
    output->Allocate();
    output->FillBuffer(0);

    ImageIteratorTypeInt classifIt(ConnectedComponentFilter->GetOutput(), ConnectedComponentFilter->GetOutput()->GetLargestPossibleRegion());
    OutputIteratorType outputIt(output, output->GetLargestPossibleRegion());
    while(!classifIt.IsAtEnd())
    {
        outputIt.Set(0);
        if(classifIt.Get()!=0)
        {
            if(enoughContour[classifIt.Get()]==true)
            {
                outputIt.Set(1);
            }
        }
        ++outputIt;
        ++classifIt;
    }
}

} //end of namespace anima
