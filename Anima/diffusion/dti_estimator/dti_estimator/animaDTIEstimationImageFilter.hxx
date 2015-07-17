#pragma once
#include "animaDTIEstimationImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>

#include <vnl/algo/vnl_determinant.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

#include <animaVectorOperations.h>

#include <boost/lexical_cast.hpp>

namespace anima
{

template <class PixelScalarType>
void
DTIEstimationImageFilter<PixelScalarType>
::AddGradientDirection(unsigned int i, vnl_vector_fixed<double,3> &grad)
{
    if (i == m_GradientDirections.size())
        m_GradientDirections.push_back(grad);
    else if (i > m_GradientDirections.size())
        std::cerr << "Trying to add a direction not contiguous... Add directions contiguously (0,1,2,3,...)..." << std::endl;
    else
        m_GradientDirections[i] = grad;
}

template <class PixelScalarType>
void
DTIEstimationImageFilter<PixelScalarType>
::CheckComputationMask()
{
    typedef itk::ImageRegionConstIterator <InputImageType> B0IteratorType;
    typedef itk::ImageRegionIterator <MaskImageType> MaskIteratorType;

    unsigned int firstB0Index = 0;
    while (m_BValuesList[firstB0Index] != 0)
        ++firstB0Index;

    B0IteratorType b0Itr(this->GetInput(firstB0Index),this->GetOutput()->GetLargestPossibleRegion());
    b0Itr.GoToBegin();

    typename MaskImageType::Pointer maskImage = MaskImageType::New();
    maskImage->Initialize();
    maskImage->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
    maskImage->SetSpacing (this->GetInput(0)->GetSpacing());
    maskImage->SetOrigin (this->GetInput(0)->GetOrigin());
    maskImage->SetDirection (this->GetInput(0)->GetDirection());
    maskImage->Allocate();

    MaskIteratorType maskItr (maskImage,this->GetOutput()->GetLargestPossibleRegion());
    while (!b0Itr.IsAtEnd())
    {
        if (b0Itr.Get() > m_B0Threshold)
            maskItr.Set(1);
        else
            maskItr.Set(0);

        ++b0Itr;
        ++maskItr;
    }

    this->SetComputationMask(maskImage);
}

template <class PixelScalarType>
void
DTIEstimationImageFilter<PixelScalarType>
::GenerateOutputInformation()
{
    // Override the method in itkImageSource, so we can set the vector length of
    // the output itk::VectorImage

    this->Superclass::GenerateOutputInformation();

    OutputImageType *output = this->GetOutput();
    output->SetVectorLength(m_NumberOfComponents);
}

template <class PixelScalarType>
void
DTIEstimationImageFilter<PixelScalarType>
::BeforeThreadedGenerateData()
{
    if (m_BValuesList.size() != this->GetNumberOfIndexedInputs())
    {
        std::string error("There should be the same number of input images and input b-values... ");
        error += boost::lexical_cast <std::string> (m_BValuesList.size());
        error += " ";
        error += boost::lexical_cast <std::string> (this->GetNumberOfIndexedInputs());

        throw itk::ExceptionObject(__FILE__, __LINE__,error,ITK_LOCATION);
    }

    itk::ImageRegionIterator <OutputImageType> fillOut(this->GetOutput(),this->GetOutput()->GetLargestPossibleRegion());
    typename OutputImageType::PixelType emptyModelVec(m_NumberOfComponents);
    emptyModelVec.Fill(0);

    while (!fillOut.IsAtEnd())
    {
        fillOut.Set(emptyModelVec);
        ++fillOut;
    }

    m_EstimatedB0Image = OutputB0ImageType::New();
    m_EstimatedB0Image->Initialize();
    m_EstimatedB0Image->SetOrigin(this->GetOutput()->GetOrigin());
    m_EstimatedB0Image->SetSpacing(this->GetOutput()->GetSpacing());
    m_EstimatedB0Image->SetDirection(this->GetOutput()->GetDirection());
    m_EstimatedB0Image->SetRegions(this->GetOutput()->GetLargestPossibleRegion());

    m_EstimatedB0Image->Allocate();
    m_EstimatedB0Image->FillBuffer(0.0);

    m_NonColinearDirectionIndexes.clear();
    for (unsigned int i = 0;i < m_BValuesList.size();++i)
    {
        if (anima::ComputeNorm(m_GradientDirections[i]) == 0)
        {
            m_NonColinearDirectionIndexes.push_back(i);
            continue;
        }

        bool colinearDataFound = false;
        for (unsigned int j = 0;j < m_NonColinearDirectionIndexes.size();++j)
        {
            if (m_BValuesList[m_NonColinearDirectionIndexes[j]] == 0)
                continue;

            if (std::abs(anima::ComputeScalarProduct(m_GradientDirections[m_NonColinearDirectionIndexes[j]],m_GradientDirections[i])) > std::cos(M_PI/180))
            {
                if (m_BValuesList[i] <= m_BValuesList[m_NonColinearDirectionIndexes[j]])
                    m_NonColinearDirectionIndexes[j] = i;

                colinearDataFound = true;
                break;
            }
        }

        if (!colinearDataFound)
            m_NonColinearDirectionIndexes.push_back(i);
    }

    m_AverageBValue = 0;
    unsigned int numBValues = 0;
    for (unsigned int i = 0;i < m_NonColinearDirectionIndexes.size();++i)
    {
        m_AverageBValue += m_BValuesList[m_NonColinearDirectionIndexes[i]];
        if (m_BValuesList[m_NonColinearDirectionIndexes[i]] != 0)
            ++numBValues;
    }

    m_AverageBValue /= numBValues;

    for (unsigned int i = 0;i < m_NonColinearDirectionIndexes.size();++i)
        m_BValuesList[m_NonColinearDirectionIndexes[i]] /= m_AverageBValue;

    // Design matrix computation

    m_DesignMatrix.set_size(m_NonColinearDirectionIndexes.size(),m_NumberOfComponents+1);

    for (unsigned int i = 0;i < m_NonColinearDirectionIndexes.size();++i)
    {
        m_DesignMatrix(i,0) = 1;
        unsigned int iIndex = m_NonColinearDirectionIndexes[i];

        unsigned int pos = 1;
        for (unsigned int j = 0;j < 3;++j)
            for (unsigned int k = 0;k <= j;++k)
            {
                if (j != k)
                    m_DesignMatrix(i,pos) = - 2 * m_BValuesList[iIndex] * m_GradientDirections[iIndex][j] * m_GradientDirections[iIndex][k];
                else
                    m_DesignMatrix(i,pos) = - m_BValuesList[iIndex] * m_GradientDirections[iIndex][j] * m_GradientDirections[iIndex][j];

                ++pos;
            }
    }

    m_SolveMatrix = vnl_matrix_inverse<double> (m_DesignMatrix.transpose() * m_DesignMatrix);
    m_SolveMatrix = m_SolveMatrix * m_DesignMatrix.transpose();

    Superclass::BeforeThreadedGenerateData();
}

template <class PixelScalarType>
void
DTIEstimationImageFilter<PixelScalarType>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionConstIterator <InputImageType> ImageIteratorType;

    unsigned int numInputs = m_NonColinearDirectionIndexes.size();
    std::vector <ImageIteratorType> inIterators;
    for (unsigned int i = 0;i < numInputs;++i)
        inIterators.push_back(ImageIteratorType(this->GetInput(m_NonColinearDirectionIndexes[i]),outputRegionForThread));

    typedef itk::ImageRegionIterator <OutputImageType> OutImageIteratorType;
    OutImageIteratorType outIterator(this->GetOutput(),outputRegionForThread);

    typedef itk::ImageRegionIteratorWithIndex <OutputB0ImageType> OutB0ImageIteratorType;
    OutB0ImageIteratorType outB0Iterator(m_EstimatedB0Image,outputRegionForThread);

    typedef typename OutputImageType::PixelType OutputPixelType;
    std::vector <double> dwi(this->GetNumberOfIndexedInputs(),0);
    OutputPixelType resVec(m_NumberOfComponents);

    vnl_matrix <double> tmpTensor(3,3);

    while (!outIterator.IsAtEnd())
    {
        for (unsigned int i = 0;i < m_NumberOfComponents;++i)
            resVec[i] = 0;

        // Compute mean B0 value and discard voxels under B0 threshold
        double b0Value = 0;
        unsigned int numberOfB0Images = 0;
        for (unsigned int j = 0;j < m_NonColinearDirectionIndexes.size();++j)
        {
            if (m_BValuesList[m_NonColinearDirectionIndexes[j]] == 0)
            {
                b0Value += inIterators[j].Get();
                ++numberOfB0Images;
            }
        }

        b0Value /= numberOfB0Images;

        if (b0Value <= m_B0Threshold)
        {
            for (unsigned int i = 0;i < numInputs;++i)
                ++inIterators[i];

            outIterator.Set(resVec);

            ++outIterator;
            ++outB0Iterator;
            continue;
        }

        for (unsigned int i = 0;i < numInputs;++i)
        {
            if (inIterators[i].Get() > 0)
                dwi[i] = log(inIterators[i].Get());
            else
                dwi[i] = log(1.0e-16);
        }

        // We don't care about the B0 log-value
        for (unsigned int i = 1;i < m_NumberOfComponents+1;++i)
        {
            for (unsigned int j = 0;j < numInputs;++j)
                resVec[i-1] += m_SolveMatrix(i,j) * dwi[j];
        }

        double outB0Value = 0;
        for (unsigned int j = 0;j < numInputs;++j)
            outB0Value += m_SolveMatrix(0,j) * dwi[j];

        outB0Value = exp(outB0Value);

        unsigned int pos = 0;
        for (unsigned int i = 0;i < 3;++i)
            for (unsigned int j = 0;j <= i;++j)
            {
                tmpTensor(i,j) = resVec[pos];
                if (i != j)
                    tmpTensor(j,i) = tmpTensor(i,j);
                ++pos;
            }

        vnl_symmetric_eigensystem <double> tmpEigs(tmpTensor);

        bool isTensorOk = true;

        if (m_RemoveDegeneratedTensors)
        {
            for (unsigned int i = 0;i < 3;++i)
            {
                if (tmpEigs.D[i] <= 0)
                {
                    isTensorOk = false;
                    break;
                }
            }
        }

        if (isTensorOk)
        {
            resVec /= m_AverageBValue;

            outIterator.Set(resVec);
            outB0Iterator.Set(outB0Value);
        }

        for (unsigned int i = 0;i < numInputs;++i)
            ++inIterators[i];

        ++outIterator;
        ++outB0Iterator;
    }
}

} // end of namespace anima
