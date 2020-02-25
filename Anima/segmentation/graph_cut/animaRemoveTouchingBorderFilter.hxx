#pragma once

#include "animaRemoveTouchingBorderFilter.h"

namespace anima
{

template<typename TInput, typename TMask, typename TOutput>
void RemoveTouchingBorderFilter< TInput, TMask, TOutput >::SetInputImageSeg(const TInput* image)
{
    this->SetNthInput(0, const_cast<TInput*>(image));
}

template<typename TInput, typename TMask, typename TOutput>
void RemoveTouchingBorderFilter< TInput, TMask, TOutput >::SetMask(const TMask* mask)
{
    this->SetNthInput(1, const_cast<TMask*>(mask));
}


template<typename TInput, typename TMask, typename TOutput>
typename TInput::ConstPointer RemoveTouchingBorderFilter< TInput, TMask, TOutput >::GetInputImageSeg()
{
    return static_cast< const TInput * >
            ( this->itk::ProcessObject::GetInput(0) );
}
template<typename TInput, typename TMask, typename TOutput>
typename TMask::ConstPointer RemoveTouchingBorderFilter< TInput, TMask, TOutput >::GetMask()
{
    return static_cast< const TMask * >
            ( this->itk::ProcessObject::GetInput(1) );
}


template<typename TInput, typename TMask, typename TOutput>
itk::DataObject::Pointer RemoveTouchingBorderFilter< TInput, TMask, TOutput >::MakeOutput(unsigned int idx)
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

template<typename TInput, typename TMask, typename TOutput>
TOutput* RemoveTouchingBorderFilter< TInput, TMask, TOutput >::GetOutputNonTouchingBorder()
{
    return dynamic_cast< TOutput* >( this->itk::ProcessObject::GetOutput(0) );
}

template<typename TInput, typename TMask, typename TOutput>
TOutput* RemoveTouchingBorderFilter< TInput, TMask, TOutput >::GetOutputTouchingBorder()
{
    return dynamic_cast< TOutput* >( this->itk::ProcessObject::GetOutput(1) );
}

template<typename TInput, typename TMask, typename TOutput>
void
RemoveTouchingBorderFilter< TInput, TMask, TOutput >::
WriteOutputs()
{
    if( m_OutputNonTouchingBorderFilename != "" )
    {
        std::cout << "Writing output non touching border image to: " << m_OutputNonTouchingBorderFilename << std::endl;
        anima::writeImage<TOutput>(m_OutputNonTouchingBorderFilename, this->GetOutputNonTouchingBorder());
    }

    if( m_OutputTouchingBorderFilename != "" )
    {
        std::cout << "Writing output touching border image to: " << m_OutputTouchingBorderFilename << std::endl;
        anima::writeImage<TOutput>(m_OutputTouchingBorderFilename, this->GetOutputTouchingBorder());
    }
}

template<typename TInput, typename TMask, typename TOutput>
void
RemoveTouchingBorderFilter< TInput, TMask, TOutput >::
GenerateData()
{
    InputImagePointer InputImageSeg = InputImageType::New();
    InputImageSeg->SetRegions(this->GetInputImageSeg()->GetLargestPossibleRegion());
    InputImageSeg->CopyInformation(this->GetInputImageSeg());
    InputImageSeg->Allocate();
    InputImageSeg->FillBuffer(0);

    // Find labels overlaping mask border
    InputIteratorType InputImageSegIt(InputImageSeg, InputImageSeg->GetLargestPossibleRegion() );
    InputConstIteratorType maskSegIt(this->GetInputImageSeg(), this->GetInputImageSeg()->GetLargestPossibleRegion() );
    
    while(!maskSegIt.IsAtEnd())
    {
        InputImageSegIt.Set(maskSegIt.Get());

        ++maskSegIt;
        ++InputImageSegIt;
    }

    // --------------------------- Lesions
    typename ConnectedComponentImageFilterType::Pointer connectedComponentImageFilter = ConnectedComponentImageFilterType::New();
    typename LabelImageToLabelMapFilterType::Pointer labelImageToLabelMapFilter = LabelImageToLabelMapFilterType::New();
    typename LabelImageToLabelMapFilterType2::Pointer labelImageToLabelMapFilter2 = LabelImageToLabelMapFilterType2::New();
    typename LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
    typename TOutputMap::Pointer labelMap = ITK_NULLPTR;
    int originalNumberOfObject = 0;
    bool connectivity = false;
    
    if( !m_LabeledImage )
    {
        // Create the labeled image of the input segmentation
        connectedComponentImageFilter->SetInput( InputImageSeg );
        connectedComponentImageFilter->SetFullyConnected( connectivity );
        connectedComponentImageFilter->SetNumberOfWorkUnits( this->GetNumberOfWorkUnits() );
        connectedComponentImageFilter->Update();

        labelImageToLabelMapFilter->SetInput( connectedComponentImageFilter->GetOutput() );
        labelImageToLabelMapFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
        labelImageToLabelMapFilter->Update(); // The output of this filter is an itk::LabelMap, which contains itk::LabelObject's
        originalNumberOfObject = labelImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() ;

        labelMap = labelImageToLabelMapFilter->GetOutput();
        labelMap->DisconnectPipeline();
    }
    else
    {
        labelImageToLabelMapFilter2->SetInput( InputImageSeg );
        labelImageToLabelMapFilter2->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
        labelImageToLabelMapFilter2->Update(); // The output of this filter is an itk::LabelMap, which contains itk::LabelObject's
        originalNumberOfObject = labelImageToLabelMapFilter2->GetOutput()->GetNumberOfLabelObjects() ;

        labelMap = labelImageToLabelMapFilter2->GetOutput();
        labelMap->DisconnectPipeline();
    }

    // Convert a LabelMap to a normal image with different values representing each region
    labelMapToLabelImageFilter->SetInput( labelMap );
    labelMapToLabelImageFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    labelMapToLabelImageFilter->Update();

    // ------------------------------------- Mask bordure

    ImageIteratorTypeInt segLabelIt (labelMapToLabelImageFilter->GetOutput(), labelMapToLabelImageFilter->GetOutput()->GetLargestPossibleRegion() );
    
    std::vector<int> labelsToRemove;
    
    if(!m_NoContour)
    {
        // Get envelop of the mask
        // Generate connected components
        typename ConnectedComponentImageFilterType2::Pointer connectedComponentImageFilter2  = ConnectedComponentImageFilterType2::New();
        connectedComponentImageFilter2->SetInput( this->GetMask() );
        connectedComponentImageFilter2->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

        // Generate contours for each component
        itk::LabelContourImageFilter<ImageTypeInt, ImageTypeInt>::Pointer labelContourImageFilter = itk::LabelContourImageFilter<ImageTypeInt, ImageTypeInt>::New();
        labelContourImageFilter->SetInput( connectedComponentImageFilter2->GetOutput() );
        labelContourImageFilter->SetFullyConnected(true);
        labelContourImageFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
        labelContourImageFilter->SetBackgroundValue(0);
        labelContourImageFilter->Update();


        // Find labels overlaping mask border
        ImageIteratorTypeInt maskContourIt (labelContourImageFilter->GetOutput(), labelContourImageFilter->GetOutput()->GetLargestPossibleRegion() );

        std::vector<int>::iterator it;
        while(!maskContourIt.IsAtEnd())
        {
            if( maskContourIt.Get()!=0 && segLabelIt.Get()!=0 )
            {
                it = find (labelsToRemove.begin(), labelsToRemove.end(), segLabelIt.Get());
                if (it == labelsToRemove.end())
                {
                    labelsToRemove.push_back( segLabelIt.Get() );
                    m_labelsToRemove.push_back( segLabelIt.Get() );
                }
            }
            ++maskContourIt;
            ++segLabelIt;
        }
    }
    else
    {
        // Find labels overlaping mask border
        MaskConstIteratorType maskContourIt (this->GetMask(), this->GetMask()->GetLargestPossibleRegion() );

        std::vector<int>::iterator it;
        while(!maskContourIt.IsAtEnd())
        {
            if( maskContourIt.Get()!=0 && segLabelIt.Get()!=0)
            {
                it = find (labelsToRemove.begin(), labelsToRemove.end(), segLabelIt.Get());
                if (it == labelsToRemove.end())
                {
                    labelsToRemove.push_back( segLabelIt.Get() );
                    m_labelsToRemove.push_back( segLabelIt.Get() );
                }
            }
            ++maskContourIt;
            ++segLabelIt;
        }
    }

    
    // Remove all regions that were marked for removal.
    for(unsigned int i = 0; i < labelsToRemove.size(); ++i)
    {
        labelMap->RemoveLabel(labelsToRemove[i]);
    }

    // Convert a LabelMap to a normal image with different values representing each region
    typename LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter2 = LabelMapToLabelImageFilterType::New();
    labelMapToLabelImageFilter2->SetInput( labelMap );
    labelMapToLabelImageFilter2->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    labelMapToLabelImageFilter2->Update();

    MaskConstIteratorType maskIt (this->GetMask(), this->GetMask()->GetLargestPossibleRegion() );
    ImageIteratorTypeInt segIt (labelMapToLabelImageFilter2->GetOutput(), labelMapToLabelImageFilter2->GetOutput()->GetLargestPossibleRegion() );
    maskSegIt.GoToBegin();

    OutputImagePointer output = this->GetOutputNonTouchingBorder();
    output->SetRegions(this->GetInputImageSeg()->GetLargestPossibleRegion());
    output->CopyInformation(this->GetInputImageSeg());
    output->Allocate();
    output->FillBuffer(0);

    OutputImagePointer output2 = this->GetOutputTouchingBorder();
    output2->SetRegions(this->GetInputImageSeg()->GetLargestPossibleRegion());
    output2->CopyInformation(this->GetInputImageSeg());
    output2->Allocate();
    output2->FillBuffer(0);

    OutputIteratorType outIt (output, output->GetLargestPossibleRegion() );
    OutputIteratorType out2It (output2, output2->GetLargestPossibleRegion() );

    while(!maskIt.IsAtEnd())
    {
        if( segIt.Get()!=0)
        {
            outIt.Set(1);
        }
        if( (segIt.Get()==0) && (maskSegIt.Get()!=0))
        {
            out2It.Set(1);
        }
        ++maskIt;
        ++segIt;
        ++outIt;
        ++out2It;
        ++maskSegIt;
    }

    if ( m_Verbose )
    {
        std::cout << " -- Rule to erase objects touching mask border: " << std::endl;
        std::cout << "    * Initial number of objects: " << originalNumberOfObject << std::endl;
        std::cout << "    * Number of rejected objects: " << labelsToRemove.size()  << std::endl;
        std::cout << "    * Number of objects after clean: " << labelMap->GetNumberOfLabelObjects() << std::endl;
        std::cout << std::endl;
    }

}


} //end of namespace anima
