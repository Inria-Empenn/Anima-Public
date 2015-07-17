#pragma once

#include <animaMaskedImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>

namespace anima
{

template <class PixelScalarType>
class DTIEstimationImageFilter :
public anima::MaskedImageToImageFilter < itk::Image<PixelScalarType,3>, itk::VectorImage<PixelScalarType,3> >
{
public:
    /** Standard class typedefs. */
    typedef DTIEstimationImageFilter<PixelScalarType> Self;
    typedef itk::Image<PixelScalarType,3> InputImageType;
    typedef itk::VectorImage<PixelScalarType,3> OutputImageType;
    typedef InputImageType OutputB0ImageType;
    typedef OutputImageType DTIImageType;
    typedef itk::Image<PixelScalarType,4> Image4DType;
    typedef anima::MaskedImageToImageFilter< InputImageType, OutputImageType > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods) */
    itkTypeMacro(DTIEstimationImageFilter, MaskedImageToImageFilter);

    /** Image typedef support */
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::MaskImageType MaskImageType;
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    void SetBValuesList(std::vector <double> bValuesList ){m_BValuesList = bValuesList;}

    itkSetMacro(B0Threshold, float);
    itkGetMacro(B0Threshold, float);

    itkSetMacro(RemoveDegeneratedTensors, bool);

    itkGetMacro(EstimatedB0Image, OutputB0ImageType *);

    void AddGradientDirection(unsigned int i, vnl_vector_fixed<double,3> &grad);

protected:
    DTIEstimationImageFilter()
        : Superclass()
    {
        m_BValuesList.clear();

        m_B0Threshold = 0;
        m_AverageBValue = 1;
        m_RemoveDegeneratedTensors = true;
        m_EstimatedB0Image = NULL;
    }

    virtual ~DTIEstimationImageFilter() {}

    void CheckComputationMask();

    void GenerateOutputInformation(void);
    void BeforeThreadedGenerateData();
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId);

    void PrintSelf(std::ostream& os, itk::Indent indent) const
    {
        Superclass::PrintSelf(os,indent);
    }

private:
    DTIEstimationImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    std::vector <double> m_BValuesList;
    std::vector< vnl_vector_fixed<double,3> > m_GradientDirections;
    std::vector <unsigned int> m_NonColinearDirectionIndexes;

    float m_B0Threshold;
    typename OutputB0ImageType::Pointer m_EstimatedB0Image;

    static const unsigned int m_NumberOfComponents = 6;
    double m_AverageBValue;

    vnl_matrix <double> m_DesignMatrix, m_SolveMatrix;
    bool m_RemoveDegeneratedTensors;
};

} // end of namespace anima

#include "animaDTIEstimationImageFilter.hxx"
