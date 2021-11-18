#pragma once

#include <animaDTIScalarMapsImageFilter.h>

#include <itkImageRegionConstIterator.h>
#include <itkNeighborhoodInnerProduct.h>
#include <itkImageRegionIterator.h>
#include <itkNeighborhoodAlgorithm.h>
#include <itkZeroFluxNeumannBoundaryCondition.h>
#include <itkOffset.h>
#include <itkProgressReporter.h>
#include <itkMatrix.h>
#include <itkSymmetricEigenAnalysis.h>

#include <animaBaseTensorTools.h>

namespace anima
{

template <unsigned int ImageDimension>
DTIScalarMapsImageFilter < ImageDimension >::DTIScalarMapsImageFilter() :
    Superclass()
{
    this->SetNumberOfRequiredOutputs(5);

    for (unsigned int i = 0;i < 5;++i)
        this->SetNthOutput(i, this->MakeOutput(i));
}

/**
 *   Make Ouput
 */
template <unsigned int ImageDimension>
itk::DataObject::Pointer
DTIScalarMapsImageFilter< ImageDimension >::MakeOutput(itk::ProcessObject::DataObjectPointerArraySizeType idx)
{
    return (OutputImageType::New()).GetPointer();
}

template <unsigned int ImageDimension>
void
DTIScalarMapsImageFilter< ImageDimension >::SetAnglesMatrix(vnl_matrix <double> &affMatrix)
{
    vnl_matrix<double> tmpMat(3,3);
    m_RigidAnglesMatrix.set_size(3,3);
    anima::ExtractRotationFromJacobianMatrix(affMatrix,m_RigidAnglesMatrix,tmpMat);
}


template <unsigned int ImageDimension>
void
DTIScalarMapsImageFilter< ImageDimension >
::DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread)
{
    itk::ImageRegionConstIterator<TensorImageType> tensorIterator;
    itk::ImageRegionIterator<OutputImageType> adcIterator;

    // Allocate output
    typename  InputImageType::ConstPointer tensorImage  = this->GetInput();
    typename OutputImageType::Pointer adcImage = this->GetOutput(0);

    TensorImageRegionType tensorRegionForThread;
    tensorRegionForThread.SetIndex(outputRegionForThread.GetIndex());
    tensorRegionForThread.SetSize(outputRegionForThread.GetSize());

    tensorIterator = itk::ImageRegionConstIterator<TensorImageType>(tensorImage, tensorRegionForThread);
    adcIterator = itk::ImageRegionIterator<OutputImageType>(adcImage, outputRegionForThread);

    itk::ImageRegionIterator<OutputImageType> faIterator;
    typename OutputImageType::Pointer faImage = this->GetOutput(1);
    faIterator = itk::ImageRegionIterator<OutputImageType>(faImage, outputRegionForThread);

    itk::ImageRegionIterator<OutputImageType> axialIterator;
    typename OutputImageType::Pointer axialImage = this->GetOutput(2);
    axialIterator = itk::ImageRegionIterator<OutputImageType>(axialImage, outputRegionForThread);

    itk::ImageRegionIterator<OutputImageType> radialIterator;
    typename OutputImageType::Pointer radialImage = this->GetOutput(3);
    radialIterator = itk::ImageRegionIterator<OutputImageType>(radialImage, outputRegionForThread);

    itk::ImageRegionIterator<OutputImageType> angleIterator;
    typename OutputImageType::Pointer angleImage = this->GetOutput(4);
    angleIterator = itk::ImageRegionIterator<OutputImageType>(angleImage, outputRegionForThread);

    typedef vnl_matrix<double> TensorSymMatrixType;
    typedef itk::Vector<double, 3> EigenVectorType;

    EigenVectorType eigenValues;
    TensorSymMatrixType tensorSymMatrix(3,3);
    TensorSymMatrixType eigenVectors(3,3);

    itk::SymmetricEigenAnalysis<TensorSymMatrixType, EigenVectorType> eigenAnalysis;
    eigenAnalysis.SetDimension(3);

    while (!tensorIterator.IsAtEnd())
    {
        TensorVectorType tensor = tensorIterator.Get();
        double a(tensor[0]), c(tensor[2]), f(tensor[5]), ADC(0);
        ADC = (a+c+f) / 3.0;

        adcIterator.Set(ADC);

        anima::GetTensorFromVectorRepresentation(tensor, tensorSymMatrix,3,false);

        eigenAnalysis.ComputeEigenValuesAndVectors(tensorSymMatrix, eigenValues, eigenVectors);

        double l1(eigenValues[2]), l2(eigenValues[1]), l3(eigenValues[0]), fa(1);
        double num = std::sqrt ((l1 -l2) * (l1 -l2) + (l2 -l3) * (l2 -l3) + (l3 - l1) * (l3 - l1));
        double den = std::sqrt (l1*l1 + l2*l2 + l3*l3);

        if (den == 0)
            fa = 0;
        else
            fa = std::sqrt(0.5) * (num / den);

        faIterator.Set(fa);

        axialIterator.Set(l1);
        radialIterator.Set((l2+l3) / 2.0);

        if (l1 > 0)
        {
            // scalar product with eigenVectors(2,:) * m_RigidAnglesMatrix * [0,0,1]'
            double cosAngleValue = 0.0;

            for (unsigned int j = 0;j < 3;++j)
                cosAngleValue += eigenVectors(2,j) * m_RigidAnglesMatrix(j,2);

            cosAngleValue = std::abs(cosAngleValue);
            angleIterator.Set(std::acos(cosAngleValue) * 180.0 / M_PI);
        }

        this->IncrementNumberOfProcessedPoints();
        ++tensorIterator;
        ++adcIterator;
        ++faIterator;
        ++axialIterator;
        ++radialIterator;
        ++angleIterator;
    }
}

} // end of namespace anima
