#pragma once

#include <animaBaseTensorTools.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

namespace anima
{

template <unsigned int ImageDimension>
void
FlipTensorsImageFilter<ImageDimension>
::GenerateOutputInformation()
{
    // Override the method in itkImageSource, so we can set the vector length of
    // the output itk::VectorImage
    
    this->Superclass::GenerateOutputInformation();
    
    OutputImageType *output = this->GetOutput();
    output->SetVectorLength(m_NumberOfComponents);
}

template <unsigned int ImageDimension>
void
FlipTensorsImageFilter<ImageDimension>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                       itk::ThreadIdType threadId)
{
    itk::ImageRegionConstIterator<InputImageType> inItr(this->GetInput(), outputRegionForThread);
    itk::ImageRegionConstIterator<MaskImageType> maskItr(this->GetComputationMask(), outputRegionForThread);
    itk::ImageRegionIterator<OutputImageType> outItr(this->GetOutput(), outputRegionForThread);

    typedef vnl_matrix<double> MatrixType;
    typedef vnl_diag_matrix<double> EigenVectorType;

    EigenVectorType eVals(3);
    MatrixType eVecs(3,3);
    MatrixType tensorSymMatrix(3,3);
    InputPixelType inTensor(6);
    OutputPixelType outTensor(6);

    itk::SymmetricEigenAnalysis<MatrixType, EigenVectorType, MatrixType> eigenAnalysis(3);

    while (!maskItr.IsAtEnd())
    {
        if (maskItr.Get() == 0)
        {
            outTensor.Fill(0.0);
            outItr.Set(outTensor);
            
            ++inItr;
            ++outItr;
            ++maskItr;
            continue;
        }
        
        inTensor = inItr.Get();
        anima::GetTensorFromVectorRepresentation(inTensor, tensorSymMatrix,3,false);

        eigenAnalysis.ComputeEigenValuesAndVectors(tensorSymMatrix, eVals, eVecs);
        
        if (m_FlipXAxis)
            for (unsigned int i = 0;i < 3;++i)
                eVecs(i,0) *= -1.0;
        
        if (m_FlipYAxis)
            for (unsigned int i = 0;i < 3;++i)
                eVecs(i,1) *= -1.0;
        
        if (m_FlipZAxis)
            for (unsigned int i = 0;i < 3;++i)
                eVecs(i,2) *= -1.0;
        
        anima::RecomposeTensor(eVals, eVecs, tensorSymMatrix);
        anima::GetVectorRepresentation(tensorSymMatrix,outTensor,m_NumberOfComponents,false);
        
        outItr.Set(outTensor);

        ++inItr;
        ++outItr;
        ++maskItr;
    }
}

} // end of namespace anima
