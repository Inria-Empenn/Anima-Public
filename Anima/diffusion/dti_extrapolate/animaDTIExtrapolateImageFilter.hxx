#pragma once
#include "animaDTIExtrapolateImageFilter.h"

#include <animaBaseTensorTools.h>
#include <itkSymmetricEigenAnalysis.h>

#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>

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
        OutputPixelType tmpVec(m_NumberOfComponents);
        vnl_matrix <double> tmpTensor(3,3);
        vnl_matrix <double> workTensor(3,3);
        typename OutputB0ImageType::PointType refPoint, tmpPoint;

        typedef itk::SymmetricEigenAnalysis < vnl_matrix <double>, vnl_diag_matrix<double>, vnl_matrix <double> > EigenAnalysisType;
        EigenAnalysisType eigenAnalysis(3);
        eigenAnalysis.SetOrderEigenValues(true);
        vnl_diag_matrix <double> eigVals(3);

        OutputImageRegionType totalRegion = this->GetOutput()->GetLargestPossibleRegion();
        OutputImageRegionType interpRegion;

        while (!inIterator.IsAtEnd())
        {
            resVec = inIterator.Get();

            anima::GetTensorFromVectorRepresentation(resVec,tmpTensor,3);
            eigenAnalysis.ComputeEigenValues(tmpTensor,eigVals);

            if (eigVals[0] > 0)
            {
                interpOutIterator.Set(resVec);
                interpOutB0Iterator.Set(b0Iterator.Get());

                ++b0Iterator;
                ++inIterator;
                ++interpOutB0Iterator;
                ++interpOutIterator;

                continue;
            }

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

                anima::GetTensorFromVectorRepresentation(internalIt.Get(),tmpTensor,3);
                eigenAnalysis.ComputeEigenValues(tmpTensor,eigVals);

                if (eigVals[0] <= 0)
                {
                    ++internalIt;
                    ++internalB0It;
                    continue;
                }

                m_EstimatedB0Image->TransformIndexToPhysicalPoint(internalIt.GetIndex(), tmpPoint);
                double pointDist = 0;
                for (unsigned int i = 0;i < 3;++i)
                    pointDist += (tmpPoint[i] - refPoint[i]) * (tmpPoint[i] - refPoint[i]);

                double weight = 1.0 / (1.0 + std::sqrt(pointDist));

                interpB0 += weight * internalB0It.Get();

                anima::GetTensorLogarithm(tmpTensor,workTensor);
                workTensor *= weight;
                anima::GetVectorRepresentation(workTensor,tmpVec,m_NumberOfComponents);

                resVec += tmpVec;
                sumWeights += weight;

                ++internalB0It;
                ++internalIt;
            }

            if (sumWeights > 0)
            {
                resVec /= sumWeights;

                anima::GetTensorFromVectorRepresentation(resVec,tmpTensor,3);
                anima::GetTensorExponential(tmpTensor,workTensor);
                anima::GetVectorRepresentation(workTensor,resVec,m_NumberOfComponents);

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
