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
	    if( maskSegIt.Get()!=0 )
	    {
		InputImageSegIt.Set(255);
	    }
	    ++maskSegIt;
	    ++InputImageSegIt;
	}

	// --------------------------- Lesions
	// Create the labeled image of the input segmentation
	typename BinaryImageToLabelMapFilterType::Pointer binaryImageToLabelMapFilter = BinaryImageToLabelMapFilterType::New();
	binaryImageToLabelMapFilter->SetInput( InputImageSeg );
	binaryImageToLabelMapFilter->SetNumberOfThreads(this->GetNumberOfThreads());
	binaryImageToLabelMapFilter->SetCoordinateTolerance(m_Tol);
	binaryImageToLabelMapFilter->SetDirectionTolerance(m_Tol);
	binaryImageToLabelMapFilter->Update(); // The output of this filter is an itk::LabelMap, which contains itk::LabelObject's
	int originalNumberOfObject = binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() ;


	// Convert a LabelMap to a normal image with different values representing each region
	typename LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
	labelMapToLabelImageFilter->SetInput( binaryImageToLabelMapFilter->GetOutput() );
	labelMapToLabelImageFilter->SetNumberOfThreads(this->GetNumberOfThreads());
	labelMapToLabelImageFilter->SetCoordinateTolerance(m_Tol);
	labelMapToLabelImageFilter->SetDirectionTolerance(m_Tol);
	labelMapToLabelImageFilter->Update();

	// ------------------------------------- Mask bordure
	// Get envelop of the mask
	// Generate connected components
	typedef itk::ConnectedComponentImageFilter <TMask, ImageTypeInt > ConnectedComponentImageFilterType;
	typename ConnectedComponentImageFilterType::Pointer connectedComponentImageFilter  = ConnectedComponentImageFilterType::New();
	connectedComponentImageFilter->SetInput( this->GetMask() );
	connectedComponentImageFilter->SetNumberOfThreads(this->GetNumberOfThreads());
	connectedComponentImageFilter->SetCoordinateTolerance(m_Tol);
	connectedComponentImageFilter->SetDirectionTolerance(m_Tol);
        
	// Generate contours for each component
	typedef itk::LabelContourImageFilter<ImageTypeInt, ImageTypeInt> LabelContourImageFilterType;
	LabelContourImageFilterType::Pointer labelContourImageFilter = LabelContourImageFilterType::New();
	labelContourImageFilter->SetInput( connectedComponentImageFilter->GetOutput() );
	labelContourImageFilter->SetFullyConnected(true);
	labelContourImageFilter->SetNumberOfThreads(this->GetNumberOfThreads());
	labelContourImageFilter->SetCoordinateTolerance(m_Tol);
	labelContourImageFilter->SetDirectionTolerance(m_Tol);
	labelContourImageFilter->Update();

	// Find labels overlaping mask border
	ImageIteratorTypeInt maskContourIt (labelContourImageFilter->GetOutput(), labelContourImageFilter->GetOutput()->GetLargestPossibleRegion() );
	    
	if( m_NoContour )
	{
	    // Find labels overlaping mask border
	    maskContourIt = ImageIteratorTypeInt (connectedComponentImageFilter->GetOutput(), connectedComponentImageFilter->GetOutput()->GetLargestPossibleRegion() );
	}
	
	ImageIteratorTypeInt segLabelIt (labelMapToLabelImageFilter->GetOutput(), labelMapToLabelImageFilter->GetOutput()->GetLargestPossibleRegion() );

	std::vector<int> labelsToRemove;
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


	// Remove all regions that were marked for removal.
	for(unsigned int i = 0; i < labelsToRemove.size(); ++i)
	{
	    binaryImageToLabelMapFilter->GetOutput()->RemoveLabel(labelsToRemove[i]);
	}

	// Convert a LabelMap to a normal image with different values representing each region
	typename LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter2 = LabelMapToLabelImageFilterType::New();
	labelMapToLabelImageFilter2->SetInput( binaryImageToLabelMapFilter->GetOutput() );
	labelMapToLabelImageFilter2->SetNumberOfThreads(this->GetNumberOfThreads());
	labelMapToLabelImageFilter2->SetCoordinateTolerance(m_Tol);
	labelMapToLabelImageFilter2->SetDirectionTolerance(m_Tol);
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
	    std::cout << "    * Intial number of objects: " << originalNumberOfObject << std::endl;
	    std::cout << "    * Number of rejected objects: " << labelsToRemove.size()  << std::endl;
	    std::cout << "    * Number of objects after clean: " << binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << std::endl;
	    std::cout << std::endl;
	}

    }
} //end of namespace anima
