#pragma once

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

#include <animaBaseTensorTools.h>

namespace anima
{
template <class TScalarType, unsigned int NDimensions>
void
ExpTensorImageFilter<TScalarType,NDimensions>::
GenerateOutputInformation()
{
    // Override the method in itkImageSource, so we can set the vector length of
    // the output itk::VectorImage

    this->Superclass::GenerateOutputInformation();

    m_VectorSize = this->GetInput(0)->GetNumberOfComponentsPerPixel();
    TOutputImage *output = this->GetOutput();
    output->SetVectorLength(m_VectorSize);
}

template <class TScalarType, unsigned int NDimensions>
void
ExpTensorImageFilter<TScalarType,NDimensions>::
BeforeThreadedGenerateData ()
{
    unsigned int nbInputs = this->GetNumberOfIndexedInputs();
    if (nbInputs != 1)
        throw itk::ExceptionObject(__FILE__, __LINE__,"There should be one input... Exiting...",ITK_LOCATION);

    m_VectorSize = this->GetInput(0)->GetNumberOfComponentsPerPixel();
    m_TensorDimension = floor( ( sqrt((float)(8 * m_VectorSize + 1)) - 1) / 2.0);

    this->GetOutput()->SetNumberOfComponentsPerPixel(m_VectorSize);
    this->AllocateOutputs();
}

template <class TScalarType, unsigned int NDimensions>
void
ExpTensorImageFilter<TScalarType,NDimensions>::
ThreadedGenerateData (const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionConstIterator< TInputImage > InIteratorType;
    typedef itk::ImageRegionIterator< TOutputImage > OutRegionIteratorType;

    InIteratorType inIterator(this->GetInput(), outputRegionForThread);
    OutRegionIteratorType outIterator(this->GetOutput(), outputRegionForThread);

    OutputPixelType outValue;
    vnl_matrix <double> tmpTensor(m_TensorDimension, m_TensorDimension);
    vnl_matrix <double> tmpExpTensor(m_TensorDimension, m_TensorDimension);

    while (!outIterator.IsAtEnd())
    {
        outValue = inIterator.Get();

        if (!isZero(outValue))
        {
            anima::GetTensorFromVectorRepresentation(outValue,tmpTensor,m_TensorDimension,m_ScaleNonDiagonal);

            anima::GetTensorExponential(tmpTensor,tmpExpTensor);

            anima::GetVectorRepresentation(tmpExpTensor,outValue,m_VectorSize);
        }

        outIterator.Set(outValue);

        ++inIterator;
        ++outIterator;
    }
}

} // end of namespace anima
