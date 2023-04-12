#pragma once

#include "animaConnectedComponentsMetricsFilter.h"

namespace anima
{
template <typename ImageType, typename OutputType>
void ConnectedComponentsMetricsFilter<ImageType, OutputType>::SetInputReference(const ImageType* image)
{
    this->SetNthInput(0, const_cast<ImageType*>(image));
}

template <typename ImageType, typename OutputType>
void ConnectedComponentsMetricsFilter<ImageType, OutputType>::SetInputTest(const ImageType* image)
{
    this->SetNthInput(1, const_cast<ImageType*>(image));
}

template <typename ImageType, typename OutputType>
typename ImageType::ConstPointer ConnectedComponentsMetricsFilter<ImageType, OutputType>::GetInputReference()
{
    return static_cast< const  ImageType * >
            ( this->itk::ProcessObject::GetInput(0) );
}

template <typename ImageType, typename OutputType>
typename ImageType::ConstPointer ConnectedComponentsMetricsFilter<ImageType, OutputType>::GetInputTest()
{
    return static_cast< const ImageType * >
            ( this->itk::ProcessObject::GetInput(1) );
}

template <typename ImageType, typename OutputType>
OutputType* ConnectedComponentsMetricsFilter<ImageType, OutputType>::GetOutputCCReference()
{
    return dynamic_cast< OutputType * >( this->itk::ProcessObject::GetOutput(0) );
}
template <typename ImageType, typename OutputType>
OutputType *ConnectedComponentsMetricsFilter<ImageType, OutputType>::GetOutputCCTest()
{
    return dynamic_cast< OutputType * >( this->itk::ProcessObject::GetOutput(1) );
}
template <typename ImageType, typename OutputType>
OutputType * ConnectedComponentsMetricsFilter<ImageType, OutputType>::GetOutputDiffVoxelWise()
{
    return dynamic_cast< OutputType * >( this->itk::ProcessObject::GetOutput(2) );
}
template <typename ImageType, typename OutputType>
OutputType *ConnectedComponentsMetricsFilter<ImageType, OutputType>::GetOutputDiffTest()
{
    return dynamic_cast< OutputType * >( this->itk::ProcessObject::GetOutput(3) );
}
template <typename ImageType, typename OutputType>
OutputType *ConnectedComponentsMetricsFilter<ImageType, OutputType>::GetOutputDiffRef()
{
    return dynamic_cast<  OutputType * >( this->itk::ProcessObject::GetOutput(4) );
}
template <typename ImageType, typename OutputType>
OutputType *ConnectedComponentsMetricsFilter<ImageType, OutputType>::GetOutputEvolutionRef()
{
    return dynamic_cast< OutputType * >( this->itk::ProcessObject::GetOutput(5) );
}
template <typename ImageType, typename OutputType>
OutputType *ConnectedComponentsMetricsFilter<ImageType, OutputType>::GetOutputEvolutionTest()
{
    return dynamic_cast< OutputType * >( this->itk::ProcessObject::GetOutput(6) );
}

template <typename ImageType, typename OutputType>
itk::DataObject::Pointer
ConnectedComponentsMetricsFilter<ImageType, OutputType>
::MakeOutput(unsigned int idx)
{
    itk::DataObject::Pointer output;
    output = (  OutputType::New() ).GetPointer();
    return output.GetPointer();
}

template <typename ImageType, typename OutputType>
void
ConnectedComponentsMetricsFilter<ImageType, OutputType>
::WriteOutputs()
{
    if( m_OutputDiffVoxelWise_filename != "" )
    {
        std::cout << "Writing output image to: " << m_OutputDiffVoxelWise_filename << std::endl;
        anima::writeImage<OutputType>(m_OutputDiffVoxelWise_filename, this->GetOutputDiffVoxelWise() );
    }
    if( m_OutputDiffTest_filename != "" )
    {
        std::cout << "Writing output image to: " << m_OutputDiffTest_filename << std::endl;
        anima::writeImage<OutputType>(m_OutputDiffTest_filename, this->GetOutputDiffTest());
    }
    if( m_OutputDiffRef_filename != "" )
    {
        std::cout << "Writing output image to: " << m_OutputDiffRef_filename << std::endl;
        anima::writeImage<OutputType>(m_OutputDiffRef_filename, this->GetOutputDiffRef());
    }
    if( m_OutputEvolutionTest_filename != "" )
    {
        std::cout << "Writing output image to: " << m_OutputEvolutionTest_filename << std::endl;
        anima::writeImage<OutputType>(m_OutputEvolutionTest_filename, this->GetOutputEvolutionTest());
    }
    if( m_OutputEvolutionRef_filename != "" )
    {
        std::cout << "Writing output image to: " << m_OutputEvolutionRef_filename << std::endl;
        anima::writeImage<OutputType>(m_OutputEvolutionRef_filename, this->GetOutputEvolutionRef());
    }
}

template <typename ImageType, typename OutputType>
void
ConnectedComponentsMetricsFilter<ImageType, OutputType>
::GenerateData()
{
    // Get information about reference image
    typename ConnectedComponentVolumeFilterType::Pointer CCReferenceFilter = ConnectedComponentVolumeFilterType::New();
    CCReferenceFilter->SetInput( this->GetInputReference() );
    CCReferenceFilter->SetFullyConnected( m_FullyConnected );
    CCReferenceFilter->SetMinSizeMM3( m_MinSizeMM3 );
    CCReferenceFilter->SetOutputFilename( m_OutputCCReference_filename );
    CCReferenceFilter->SetOutputVolumeFilename( m_volume_ref_text_filename );
    CCReferenceFilter->SetDimension( m_Dimension );
    CCReferenceFilter->SetVerbose( false );
    if(this->GetNumberOfWorkUnits() > 0)
        CCReferenceFilter->SetNumberOfWorkUnits( this->GetNumberOfWorkUnits() );
    CCReferenceFilter->GraftOutput( this->GetOutputCCReference() );
    CCReferenceFilter->Update();
    this->GraftNthOutput( 0 , CCReferenceFilter->GetOutput() );
    CCReferenceFilter->WriteOutputs();

    m_OriginalNumberOfObjectsRef = CCReferenceFilter->GetOriginalNumberOfObjects();
    m_NumberOfObjectsRef = CCReferenceFilter->GetNumberOfObjects();
    m_ReferenceTotalVolume = CCReferenceFilter->GetTotalVolume();
    m_vect_volume_ref = CCReferenceFilter->GetVolumeVector();

    // Get information about test image
    typename ConnectedComponentVolumeFilterType::Pointer CCTestFilter = ConnectedComponentVolumeFilterType::New();
    CCTestFilter->SetInput( this->GetInputTest() );
    CCTestFilter->SetFullyConnected( m_FullyConnected );
    CCTestFilter->SetMinSizeMM3( m_MinSizeMM3 );
    CCTestFilter->SetOutputFilename( m_OutputCCTest_filename );
    CCTestFilter->SetOutputVolumeFilename( m_volume_test_text_filename );
    CCTestFilter->SetDimension( m_Dimension );
    CCTestFilter->SetVerbose( false );
    if(this->GetNumberOfWorkUnits() > 0)
        CCTestFilter->SetNumberOfWorkUnits( this->GetNumberOfWorkUnits() );
    CCTestFilter->GraftOutput( this->GetOutputCCTest() );
    CCTestFilter->Update();
    this->GraftNthOutput( 1 , CCTestFilter->GetOutput() );
    CCTestFilter->WriteOutputs();

    m_OriginalNumberOfObjectsTest = CCTestFilter->GetOriginalNumberOfObjects();
    m_NumberOfObjectsTest = CCTestFilter->GetNumberOfObjects();
    m_TestTotalVolume = CCTestFilter->GetTotalVolume();
    m_vect_volume_test = CCTestFilter->GetVolumeVector();

    this->ComputeIntersectionBetweenCC();
    this->ComputeDetection();
    this->ComputeEvolution();
    this->ComputeMetrics();
    this->PrintImagesOutput();

    if(m_metrics_text_filename.size()!=0)
    {
        std::ofstream fp(m_metrics_text_filename.c_str(), std::ios_base::out | std::ios_base::trunc);
        if(!fp)
        {
            std::cerr << "cannot open output file" << m_metrics_text_filename << std::endl;
        }
        else
        {
            this->Print(fp);
            fp.close();
        }
    }
    if(m_Verbose)
    {
        this->Print(std::cout);
    }
}

template <typename ImageType, typename OutputType>
void
ConnectedComponentsMetricsFilter<ImageType, OutputType>
::ComputeIntersectionBetweenCC()
{
    m_IntersectionVoxels.resize(m_NumberOfObjectsRef+1, std::vector<unsigned short int> (m_NumberOfObjectsTest+1,0));
    m_IntersectionMM3.resize(m_NumberOfObjectsRef+1, std::vector<double> (m_NumberOfObjectsTest+1,0));
    m_RefInclusDansTest.resize(m_NumberOfObjectsRef+1, std::vector<unsigned short int> (m_NumberOfObjectsTest+1,0));
    m_TestInclusDansRef.resize(m_NumberOfObjectsRef+1, std::vector<unsigned short int> (m_NumberOfObjectsTest+1,0));

    // Compute number of intersection voxels between each connected components
    OutputIteratorType ccReferenceImageIt(this->GetOutputCCReference(), this->GetOutputCCReference()->GetLargestPossibleRegion());
    OutputIteratorType ccTestImageIt(this->GetOutputCCTest(), this->GetOutputCCTest()->GetLargestPossibleRegion());
    while(!ccReferenceImageIt.IsAtEnd())
    {
        if((ccReferenceImageIt.Get()!=0) && (ccTestImageIt.Get()!=0))
        {
            (m_IntersectionVoxels[ccReferenceImageIt.Get()][ccTestImageIt.Get()])++;
        }
        ++ccReferenceImageIt;
        ++ccTestImageIt;
    }

    // Compute intersection between each connected components in mm3
    // Identify which connected components intersect enough (according to alpha ratio) with others.
    for (unsigned int i = 1; i < m_NumberOfObjectsRef+1; ++i)
    {
        for (unsigned int j = 1; j < m_NumberOfObjectsTest+1; ++j)
        {
            m_IntersectionMM3[i][j] = static_cast<double>(m_IntersectionVoxels[i][j]) * m_SpacingTot;
            if( ( m_IntersectionMM3[i][j] / m_vect_volume_test[j] ) > m_Alpha)
            {
                m_TestInclusDansRef[i][j] = 1;
            }
            if( ( m_IntersectionMM3[i][j] / m_vect_volume_ref[i] ) > m_Alpha)
            {
                m_RefInclusDansTest[i][j] = 1;
            }
        }
    }
}

template <typename ImageType, typename OutputType>
void
ConnectedComponentsMetricsFilter<ImageType, OutputType>
::ComputeDetection()
{
    m_vect_true_positive_ref.resize(m_NumberOfObjectsRef+1,0);
    m_vect_true_positive_test.resize(m_NumberOfObjectsTest+1,0);
    m_vect_false_positive_test.resize(m_NumberOfObjectsTest+1,0);

    for(unsigned int i = 1; i < m_NumberOfObjectsRef+1; i++)
    {
        for(unsigned int j = 1; j < m_NumberOfObjectsTest+1; j++)
        {
            if(m_OverlapDetectionType)
            {
                if(m_TestInclusDansRef[i][j] || m_RefInclusDansTest[i][j])
                {
                    m_vect_true_positive_ref[i] = 1;
                    m_vect_true_positive_test[j] = 1;
                }
            }
            else
            {
                if(m_IntersectionVoxels[i][j])
                {
                    m_vect_true_positive_ref[i] = 1;
                    m_vect_true_positive_test[j] = 1;
                }
            }
        }
    }

    m_nb_true_positive_test = 0;
    m_nb_true_positive_ref = 0;
    m_nb_false_positive = 0;
    m_nb_false_negative = 0;

    for(unsigned int i = 1; i < m_NumberOfObjectsRef+1; i++)
    {
        if(m_vect_true_positive_ref[i])
        {
            m_nb_true_positive_ref++;
        }
        else
        {
            m_nb_false_negative++;
        }
    }

    for(unsigned int j = 1; j < m_NumberOfObjectsTest+1; j++)
    {
        if(m_vect_true_positive_test[j])
        {
            m_nb_true_positive_test++;
        }
        else
        {
            m_nb_false_positive++;
            m_vect_false_positive_test[j] = 1;
        }
    }
}


template <typename ImageType, typename OutputType>
void
ConnectedComponentsMetricsFilter<ImageType, OutputType>
::ComputeEvolution()
{
    std::vector<double> ref_asso_vol(m_NumberOfObjectsRef+1,0);
    m_vect_shrink_ref.resize(m_NumberOfObjectsRef+1,0);
    m_vect_grow_ref.resize(m_NumberOfObjectsRef+1,0);
    m_vect_same_ref.resize(m_NumberOfObjectsRef+1,0);

    for(unsigned int i = 1; i < m_NumberOfObjectsRef+1; i++)
    {
        for(unsigned int j = 1; j < m_NumberOfObjectsTest+1; j++)
        {
            if(m_OverlapDetectionType)
            {
                if(m_TestInclusDansRef[i][j] || m_RefInclusDansTest[i][j])
                {
                    ref_asso_vol[i] += m_vect_volume_test[j];
                }
            }
            else
            {
                if(m_IntersectionVoxels[i][j])
                {
                    ref_asso_vol[i] += m_vect_volume_test[j];
                }
            }
        }

        if(ref_asso_vol[i]!=0)
        {
            if(m_vect_volume_ref[i]*(1+m_Beta) <= ref_asso_vol[i])
            {
                m_vect_grow_ref[i] = 1;
            }
            else if(m_vect_volume_ref[i]*(1-m_Beta) >= ref_asso_vol[i])
            {
                m_vect_shrink_ref[i] = 1;
            }
            else
            {
                m_vect_same_ref[i] = 1;
            }
        }
    }

    for(unsigned int i = 1; i < m_NumberOfObjectsRef+1; i++)
    {
        if(m_vect_grow_ref[i])
        {
            m_nb_grow_ref++;
        }
        if(m_vect_shrink_ref[i])
        {
            m_nb_shrink_ref++;
        }
        if(m_vect_same_ref[i])
        {
            m_nb_same_ref++;
        }
    }

    std::vector<double> test_asso_vol(m_NumberOfObjectsTest+1,0);
    m_vect_shrink_test.resize(m_NumberOfObjectsTest+1,0);
    m_vect_grow_test.resize(m_NumberOfObjectsTest+1,0);
    m_vect_same_test.resize(m_NumberOfObjectsTest+1,0);

    for(unsigned int j = 1; j < m_NumberOfObjectsTest+1; j++)
    {
        for(unsigned int i = 1; i < m_NumberOfObjectsRef+1; i++)
        {
            if(m_OverlapDetectionType)
            {
                if(m_TestInclusDansRef[i][j] || m_RefInclusDansTest[i][j])
                {
                    test_asso_vol[j] += m_vect_volume_ref[i];
                }

            }
            else
            {
                if(m_IntersectionVoxels[i][j])
                {
                    test_asso_vol[j] += m_vect_volume_ref[i];
                }
            }
        }

        if(test_asso_vol[j]!=0)
        {
            if(m_vect_volume_test[j] >= test_asso_vol[j]*(1+m_Beta))
            {
                m_vect_grow_test[j] = 1;
            }
            else if(m_vect_volume_test[j] <= test_asso_vol[j]*(1-m_Beta))
            {
                m_vect_shrink_test[j] = 1;
            }
            else
            {
                m_vect_same_test[j] = 1;
            }
        }
    }


    for(unsigned int j = 1; j < m_NumberOfObjectsTest+1; j++)
    {
        if(m_vect_grow_test[j])
        {
            m_nb_grow_test++;
        }
        if(m_vect_shrink_test[j])
        {
            m_nb_shrink_test++;
        }
        if(m_vect_same_test[j])
        {
            m_nb_same_test++;
        }
    }
}

template <typename ImageType, typename OutputType>
void
ConnectedComponentsMetricsFilter<ImageType, OutputType>
::ComputeMetrics()
{
    // ---------- Volume metrics ---------- //
    m_VolumeDifference = std::abs(m_ReferenceTotalVolume - m_TestTotalVolume);
    if(m_ReferenceTotalVolume!=0)
    {
        m_VolumeDifferenceRatio = m_VolumeDifference / m_ReferenceTotalVolume;
    }

    // ---------- Object detection metrics ---------- //
    if(m_NumberOfObjectsTest!=0)
    {
        m_FalsePositiveRate = static_cast<double>(m_nb_false_positive) / static_cast<double>(m_NumberOfObjectsTest);
        m_TestTruePositiveRate = static_cast<double>(m_nb_true_positive_test) / static_cast<double>(m_NumberOfObjectsTest);
    }
    if(m_NumberOfObjectsRef!=0)
    {
        m_FalseNegativeRate = static_cast<double>(m_nb_false_negative) / static_cast<double>(m_NumberOfObjectsRef);
        m_RefTruePositiveRate = static_cast<double>(m_nb_true_positive_ref) / static_cast<double>(m_NumberOfObjectsRef);
    }

    // ---------- Object evolution metrics ---------- //
    if(m_NumberOfObjectsTest!=0)
    {
        m_TestGrowRate = static_cast<double>(m_nb_grow_test) / static_cast<double>(m_NumberOfObjectsTest);
        m_TestShrinkRate = static_cast<double>(m_nb_shrink_test) / static_cast<double>(m_NumberOfObjectsTest);
        m_TestSameRate = static_cast<double>(m_nb_same_test) / static_cast<double>(m_NumberOfObjectsTest);
    }
    if(m_NumberOfObjectsRef!=0)
    {
        m_RefGrowRate = static_cast<double>(m_nb_grow_ref) / static_cast<double>(m_NumberOfObjectsRef);
        m_RefShrinkRate = static_cast<double>(m_nb_shrink_ref) / static_cast<double>(m_NumberOfObjectsRef);
        m_RefSameRate = static_cast<double>(m_nb_same_ref) / static_cast<double>(m_NumberOfObjectsRef);
    }
}

template <typename ImageType, typename OutputType>
void
ConnectedComponentsMetricsFilter<ImageType, OutputType>
::PrintImagesOutput()
{
    OutputImagePointerType outputDiffVoxelWise = this->GetOutputDiffVoxelWise();
    outputDiffVoxelWise->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
    outputDiffVoxelWise->CopyInformation(this->GetInput(0));
    outputDiffVoxelWise->Allocate();
    outputDiffVoxelWise->FillBuffer(0);

    OutputImagePointerType outputDiffTest = this->GetOutputDiffTest();
    outputDiffTest->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
    outputDiffTest->CopyInformation(this->GetInput(0));
    outputDiffTest->Allocate();
    outputDiffTest->FillBuffer(0);

    OutputImagePointerType outputDiffReference = this->GetOutputDiffRef();
    outputDiffReference->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
    outputDiffReference->CopyInformation(this->GetInput(0));
    outputDiffReference->Allocate();
    outputDiffReference->FillBuffer(0);

    OutputImagePointerType outputEvolutionTest = this->GetOutputEvolutionTest();
    outputEvolutionTest->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
    outputEvolutionTest->CopyInformation(this->GetInput(0));
    outputEvolutionTest->Allocate();
    outputEvolutionTest->FillBuffer(0);

    OutputImagePointerType outputEvolutionReference = this->GetOutputEvolutionRef();
    outputEvolutionReference->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
    outputEvolutionReference->CopyInformation(this->GetInput(0));
    outputEvolutionReference->Allocate();
    outputEvolutionReference->FillBuffer(0);

    OutputIteratorType ccReferenceImageIt(this->GetOutputCCReference(), this->GetOutputCCReference()->GetLargestPossibleRegion());
    OutputIteratorType ccTestImageIt(this->GetOutputCCTest(), this->GetOutputCCTest()->GetLargestPossibleRegion());

    OutputIteratorType outputVoxelIt(outputDiffVoxelWise, outputDiffVoxelWise->GetLargestPossibleRegion());

    OutputIteratorType outputTestIt(outputDiffTest, outputDiffTest->GetLargestPossibleRegion());
    OutputIteratorType outputRefIt(outputDiffReference, outputDiffReference->GetLargestPossibleRegion());

    OutputIteratorType outputEvolTestIt(outputEvolutionTest, outputEvolutionTest->GetLargestPossibleRegion());
    OutputIteratorType outputEvolRefIt(outputEvolutionReference, outputEvolutionReference->GetLargestPossibleRegion());

    while(!ccReferenceImageIt.IsAtEnd())
    {
        if((ccReferenceImageIt.Get()!=0) && (ccTestImageIt.Get()!=0))
        {
            outputVoxelIt.Set(0); // TP
        }
        if((ccReferenceImageIt.Get()==0) && (ccTestImageIt.Get()!=0))
        {
            outputVoxelIt.Set(1); // FP
        }
        if((ccReferenceImageIt.Get()!=0) && (ccTestImageIt.Get()==0))
        {
            outputVoxelIt.Set(3); // FN
        }

        if(m_vect_false_positive_test[ccTestImageIt.Get()])
        {
            outputVoxelIt.Set(2);
        }

        if(ccReferenceImageIt.Get()!=0)
        {
            if(m_vect_true_positive_ref[ccReferenceImageIt.Get()])
            {
                outputRefIt.Set(1); // TP
                if(m_vect_grow_ref[ccReferenceImageIt.Get()])
                {
                    outputEvolRefIt.Set(4);
                }
                if(m_vect_shrink_ref[ccReferenceImageIt.Get()])
                {
                    outputEvolRefIt.Set(3);
                }
                if(m_vect_same_ref[ccReferenceImageIt.Get()])
                {
                    outputEvolRefIt.Set(1);
                }
            }
            else
            {
                outputRefIt.Set(2); // FN
                outputEvolRefIt.Set(2);
            }
        }

        if(ccTestImageIt.Get()!=0)
        {
            if(m_vect_true_positive_test[ccTestImageIt.Get()])
            {
                outputTestIt.Set(1); // TP
                if(m_vect_grow_test[ccTestImageIt.Get()])
                {
                    outputEvolTestIt.Set(4);
                }
                if(m_vect_shrink_test[ccTestImageIt.Get()])
                {
                    outputEvolTestIt.Set(3);
                }
                if(m_vect_same_test[ccTestImageIt.Get()])
                {
                    outputEvolTestIt.Set(1);
                }
            }
            else
            {
                outputTestIt.Set(2); // FP
                outputEvolTestIt.Set(2);
            }
        }

        ++ccReferenceImageIt;
        ++ccTestImageIt;
        ++outputVoxelIt;
        ++outputEvolRefIt;
        ++outputEvolTestIt;
        ++outputTestIt;
        ++outputRefIt;
    }
}



template <typename ImageType, typename OutputType>
void
ConnectedComponentsMetricsFilter<ImageType, OutputType>
::Print(std::ostream& fp)
{

    fp << "** VOLUME INFORMATION: " << std::endl;
    fp << "* In Reference Image: " << std::endl;
    fp << "-- Original number of connected components: " << m_OriginalNumberOfObjectsRef << std::endl;
    fp << "-- Number of connected components after removing small ones: " << m_NumberOfObjectsRef  << std::endl;
    fp << "-- Total volume: " << m_ReferenceTotalVolume << " mm3" << std::endl;
    fp << std::endl;

    fp << "* In Test Image: " << std::endl;
    fp << "-- Original number of connected components: " << m_OriginalNumberOfObjectsTest << std::endl;
    fp << "-- Number of connected components after removing small ones: " << m_NumberOfObjectsTest  << std::endl;
    fp << "-- Total volume: " << m_TestTotalVolume << " mm3" << std::endl;
    fp << std::endl;

    fp << "-- Volume difference: " << m_VolumeDifference << " mm3" << std::endl;
    fp << "-- Volume difference ratio: " << m_VolumeDifferenceRatio << std::endl;
    fp << std::endl;


    fp << "** DETECTION INFORMATION" << std::endl;

    fp << "* In Reference Image: " << std::endl;
    fp << "-- Number of true positives: " << m_nb_true_positive_ref << std::endl;
    fp << " (connected components number: ";
    for(unsigned int k = 1; k < m_vect_true_positive_ref.size(); k++)
    {
        if(m_vect_true_positive_ref[k])
            fp << k << " ";
    }
    fp << ")"<< std::endl;
    fp << "-- True positive rate: " << m_RefTruePositiveRate << std::endl;
    fp << std::endl;

    fp << "-- Number of false negatives (disappeared lesions): " << m_nb_false_negative << std::endl;
    fp << " (connected components number: ";
    for(unsigned int k = 1; k < m_vect_true_positive_ref.size(); k++)
    {
        if(!m_vect_true_positive_ref[k])
            fp << k << " ";
    }
    fp << ")"<< std::endl;
    fp << "-- False negatives rate: " << m_FalseNegativeRate << std::endl;
    fp << std::endl;


    fp << "* In Test Image: " << std::endl;
    fp << "-- Number of true positives: " << m_nb_true_positive_test << std::endl;
    fp << " (connected components number: ";
    for(unsigned int k = 1; k < m_vect_true_positive_test.size(); k++)
    {
        if(m_vect_true_positive_test[k])
            fp << k << " ";
    }
    fp << ")"<< std::endl;
    fp << "-- True positive rate: " << m_TestTruePositiveRate << std::endl;
    fp << std::endl;

    fp << "-- Number of false positives (new lesions): " << m_nb_false_positive << std::endl;
    fp << " (connected components number: ";
    for(unsigned int k = 1; k < m_vect_true_positive_test.size(); k++)
    {
        if(!m_vect_true_positive_test[k])
            fp << k << " ";
    }
    fp << ")"<< std::endl;
    fp << "-- False positive rate: " << m_FalsePositiveRate << std::endl;
    fp << std::endl;


    fp << "** EVOLUTION INFORMATION" << std::endl;

    fp << "* In Reference Image: " << std::endl;
    fp << "-- Number of connected components that grow: " << m_nb_grow_ref << std::endl;
    fp << " (connected components number: ";
    for(unsigned int j = 1; j < m_vect_grow_ref.size(); j++)
    {
        if(m_vect_grow_ref[j]){fp << j << " ";}
    }
    fp << ")"<< std::endl;
    fp << "-- Growing connected component rate: " << m_RefGrowRate << std::endl;
    fp << std::endl;
    fp << "-- Number of connected components that shrink: " << m_nb_shrink_ref << std::endl;
    fp << " (connected components number: ";
    for(unsigned int j = 1; j < m_vect_shrink_ref.size(); j++)
    {
        if(m_vect_shrink_ref[j]){fp << j << " ";}
    }
    fp << ")"<< std::endl;
    fp << "-- Shrinking connected component rate: " << m_RefShrinkRate << std::endl;
    fp << std::endl;
    fp << "-- Number of connected components that stay the same: " << m_nb_same_ref << std::endl;
    fp << " (connected components number: ";
    for(unsigned int j = 1; j < m_vect_same_ref.size(); j++)
    {
        if(m_vect_same_ref[j]){fp << j << " ";}
    }
    fp << ")"<< std::endl;
    fp << "-- Same connected component rate: " << m_RefSameRate << std::endl;
    fp << std::endl;


    fp << "* In Test Image: " << std::endl;
    fp << "-- Number of connected components that have grown: " << m_nb_grow_test << std::endl;
    fp << " (connected components number: ";
    for(unsigned int j = 1; j < m_vect_grow_test.size(); j++)
    {
        if(m_vect_grow_test[j]){fp << j << " ";}
    }
    fp << ")"<< std::endl;
    fp << "-- Growing connected component rate: " << m_TestGrowRate << std::endl;
    fp << std::endl;
    fp << "-- Number of connected components that have shrink: " << m_nb_shrink_test << std::endl;
    fp << " (connected components number: ";
    for(unsigned int j = 1; j < m_vect_shrink_test.size(); j++)
    {
        if(m_vect_shrink_test[j]){fp << j << " ";}
    }
    fp << ")"<< std::endl;
    fp << "-- Shrinking connected component rate: " << m_TestShrinkRate << std::endl;
    fp << std::endl;
    fp << "-- Number of connected components that have stayed the same: " << m_nb_same_test << std::endl;
    fp << " (connected components number: ";
    for(unsigned int j = 1; j < m_vect_same_test.size(); j++)
    {
        if(m_vect_same_test[j]){fp << j << " ";}
    }
    fp << ")"<< std::endl;
    fp << "-- Same connected component rate: " << m_TestSameRate << std::endl;
    fp << std::endl;
    fp << std::endl;
    fp << std::endl;

}

} //end of namespace anima
