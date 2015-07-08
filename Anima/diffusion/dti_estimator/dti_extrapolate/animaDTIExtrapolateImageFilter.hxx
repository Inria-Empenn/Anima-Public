#pragma once
#include "animaDTIExtrapolateImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>

#include <vnl/algo/vnl_determinant.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

namespace anima
{

    template <class PixelScalarType>
    void
    DTIExtrapolateImageFilter<PixelScalarType>
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
    DTIExtrapolateImageFilter<PixelScalarType>
    ::BeforeThreadedGenerateData()
    {
        Superclass::BeforeThreadedGenerateData();

        if (m_InitialEstimatedB0Image.IsNull())
            throw itk::ExceptionObject(__FILE__, __LINE__,"No initial B0 image provided",ITK_LOCATION);

        if (this->GetInput()->GetNumberOfComponentsPerPixel() != m_NumberOfComponents)
            throw itk::ExceptionObject(__FILE__, __LINE__,"Incorrect initial DTI image provided",ITK_LOCATION);

        m_EstimatedB0Image = OutputB0ImageType::New();
        m_EstimatedB0Image->Initialize();
        m_EstimatedB0Image->SetOrigin(this->GetOutput()->GetOrigin());
        m_EstimatedB0Image->SetSpacing(this->GetOutput()->GetSpacing());
        m_EstimatedB0Image->SetDirection(this->GetOutput()->GetDirection());
        m_EstimatedB0Image->SetRegions(this->GetOutput()->GetLargestPossibleRegion());

        m_EstimatedB0Image->Allocate();
        m_EstimatedB0Image->FillBuffer(0.0);
    }

    template <class PixelScalarType>
    void
    DTIExtrapolateImageFilter<PixelScalarType>
    ::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
    {
        typedef itk::ImageRegionConstIteratorWithIndex <OutputImageType> InputImageIteratorType;
        InputImageIteratorType inIterator(this->GetInput(),outputRegionForThread);

        typedef itk::ImageRegionIteratorWithIndex <OutputImageType> OutImageIteratorType;
        OutImageIteratorType interpOutIterator(this->GetOutput(),outputRegionForThread);

        typedef itk::ImageRegionIteratorWithIndex <OutputB0ImageType> OutB0ImageIteratorType;
        OutB0ImageIteratorType b0Iterator(m_InitialEstimatedB0Image,outputRegionForThread);
        OutB0ImageIteratorType interpOutB0Iterator(m_EstimatedB0Image,outputRegionForThread);

        typedef typename OutputImageType::PixelType OutputPixelType;
        OutputPixelType resVec(m_NumberOfComponents);
        vnl_matrix <double> tmpTensor(3,3);
        typename OutputB0ImageType::PointType refPoint, tmpPoint;

        OutputImageRegionType totalRegion = this->GetOutput()->GetLargestPossibleRegion();

        while (!inIterator.IsAtEnd())
        {
            resVec = inIterator.Get();

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

            for (unsigned int i = 0;i < 3;++i)
            {
                if (tmpEigs.D[i] <= 0)
                {
                    isTensorOk = false;
                    break;
                }
            }

            if (isTensorOk)
            {
                interpOutIterator.Set(resVec);
                interpOutB0Iterator.Set(b0Iterator.Get());

                ++b0Iterator;
                ++inIterator;
                ++interpOutB0Iterator;
                ++interpOutIterator;

                continue;
            }

            OutputImageRegionType interpRegion;
            for (unsigned int i = 0;i < 3;++i)
            {
                interpRegion.SetIndex(i,std::max((int)0, (int)inIterator.GetIndex()[i] - 1));
                unsigned int iFinalIndexA = inIterator.GetIndex()[i] + 1;
                unsigned int iFinalIndexB = totalRegion.GetIndex()[i] + totalRegion.GetSize()[i] - 1;
                unsigned int tmpSize = std::min(iFinalIndexA, iFinalIndexB) - interpRegion.GetIndex()[i] + 1;
                interpRegion.SetSize(i,tmpSize);
            }

            InputImageIteratorType internalIt(this->GetInput(), interpRegion);
            OutB0ImageIteratorType internalB0It(m_InitialEstimatedB0Image, interpRegion);

            resVec.Fill(0);
            double sumWeights = 0;
            double interpB0 = 0;
            m_EstimatedB0Image->TransformIndexToPhysicalPoint(inIterator.GetIndex(), refPoint);

            while (!internalIt.IsAtEnd())
            {
                if (internalIt.GetIndex() == inIterator.GetIndex())
                {
                    ++internalIt;
                    ++internalB0It;
                    continue;
                }

                pos = 0;
                for (unsigned int i = 0;i < 3;++i)
                    for (unsigned int j = 0;j <= i;++j)
                    {
                        tmpTensor(i,j) = internalIt.Get()[pos];
                        if (i != j)
                            tmpTensor(j,i) = tmpTensor(i,j);
                        ++pos;
                    }

                tmpEigs = vnl_symmetric_eigensystem <double> (tmpTensor);
                isTensorOk = true;

                for (unsigned int i = 0;i < 3;++i)
                {
                    if (tmpEigs.D[i] <= 0)
                    {
                        isTensorOk = false;
                        break;
                    }
                }

                if (!isTensorOk)
                {
                    ++internalIt;
                    ++internalB0It;
                    continue;
                }

                m_EstimatedB0Image->TransformIndexToPhysicalPoint(internalIt.GetIndex(), tmpPoint);
                double pointDist = 0;
                for (unsigned int i = 0;i < 3;++i)
                    pointDist += (tmpPoint[i] - refPoint[i]) * (tmpPoint[i] - refPoint[i]);

                double weight = 1.0 / (1.0 + sqrt(pointDist));

                interpB0 += weight * internalB0It.Get();

                for (unsigned int i = 0;i < 3;++i)
                    tmpEigs.D[i] = std::log(tmpEigs.D[i]);

                tmpTensor = tmpEigs.recompose();
                pos = 0;
                for (unsigned int i = 0;i < 3;++i)
                    for (unsigned int j = 0;j <= i;++j)
                    {
                        resVec[pos] += tmpTensor(i,j) * weight;
                        ++pos;
                    }

                sumWeights += weight;

                ++internalB0It;
                ++internalIt;
            }

            if (sumWeights > 0)
            {
                resVec /= sumWeights;
                pos = 0;
                for (unsigned int i = 0;i < 3;++i)
                    for (unsigned int j = 0;j <= i;++j)
                    {
                        tmpTensor(i,j) = resVec[pos];
                        if (i != j)
                            tmpTensor(j,i) = tmpTensor(i,j);
                        ++pos;
                    }

                tmpEigs = vnl_symmetric_eigensystem <double> (tmpTensor);

                for (unsigned int i = 0;i < 3;++i)
                    tmpEigs.D[i] = std::exp(tmpEigs.D[i]);

                tmpTensor = tmpEigs.recompose();
                pos = 0;
                for (unsigned int i = 0;i < 3;++i)
                    for (unsigned int j = 0;j <= i;++j)
                    {
                        resVec[pos] = tmpTensor(i,j);
                        ++pos;
                    }

                interpOutIterator.Set(resVec);
                interpOutB0Iterator.Set(interpB0 / sumWeights);
            }
            else
            {
                interpOutIterator.Set(inIterator.Get());
                interpOutB0Iterator.Set(b0Iterator.Get());
            }

            ++interpOutIterator;
            ++interpOutB0Iterator;
            ++inIterator;
            ++b0Iterator;
        }
    }

} // end namespace anima
