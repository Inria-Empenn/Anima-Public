#pragma once

#include <iostream>
#include <animaMaskedImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>

#include <vector>

namespace anima
{

template <class PixelScalarType>
class LocalPatchCovarianceDistanceImageFilter :
public anima::MaskedImageToImageFilter< itk::VectorImage <PixelScalarType, 3> , itk::Image <PixelScalarType, 3> >
{
public:
    /** Standard class typedefs. */
    typedef LocalPatchCovarianceDistanceImageFilter<PixelScalarType> Self;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods) */
    itkTypeMacro(LocalPatchCovarianceDistanceImageFilter, MaskedImageToImageFilter);

    /** Image typedef support */
    typedef itk::VectorImage <PixelScalarType, 3> InputImageType;
    typedef itk::Image <PixelScalarType, 3> OutputImageType;

    typedef typename InputImageType::PixelType VectorType;

    typedef vnl_matrix <double> CovarianceType;

    typedef anima::MaskedImageToImageFilter< InputImageType, OutputImageType > Superclass;
    typedef typename Superclass::MaskImageType MaskImageType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename InputImageType::IndexType InputImageIndexType;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
    itkSetMacro(PatchHalfSize, unsigned int);

protected:
    LocalPatchCovarianceDistanceImageFilter()
    : Superclass()
    {
        this->SetNumberOfRequiredOutputs(2);
        this->SetNthOutput(0,this->MakeOutput(0));
        this->SetNthOutput(1,this->MakeOutput(1));

        m_PatchHalfSize = 1;
    }

    virtual ~LocalPatchCovarianceDistanceImageFilter() {}

    void BeforeThreadedGenerateData(void);
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId);

    void PrintSelf(std::ostream& os, itk::Indent indent) const
    {
        Superclass::PrintSelf(os,indent);
    }

private:
    LocalPatchCovarianceDistanceImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    unsigned int m_PatchHalfSize;
};

} // end namespace anima

#include "animaLocalPatchCovarianceDistanceImageFilter.hxx"
