#pragma once

#include <animaMaskedImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>

#include <animaTensorCompartment.h>

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
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(DTIEstimationImageFilter, MaskedImageToImageFilter)

    /** Image typedef support */
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::MaskImageType MaskImageType;
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    struct OptimizationDataStructure
    {
        Self *filter;
        std::vector <double> dwi, predictedValues;
        anima::TensorCompartment::Pointer tensorCompartment;
    };

    void SetBValuesList(std::vector <double> bValuesList ) {m_BValuesList = bValuesList;}

    itkSetMacro(B0Threshold, double)
    itkGetMacro(B0Threshold, double)

    itkGetMacro(EstimatedB0Image, OutputB0ImageType *)

    void AddGradientDirection(unsigned int i, vnl_vector_fixed<double,3> &grad);

protected:
    DTIEstimationImageFilter()
        : Superclass()
    {
        m_BValuesList.clear();

        m_B0Threshold = 0;
        m_EstimatedB0Image = NULL;
    }

    virtual ~DTIEstimationImageFilter() {}

    void CheckComputationMask() ITK_OVERRIDE;

    void GenerateOutputInformation() ITK_OVERRIDE;
    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId) ITK_OVERRIDE;

    static double OptimizationFunction(const std::vector<double> &x, std::vector<double> &grad, void *func_data);
    double ComputeCostAtPosition(const std::vector<double> &x, std::vector<double> &funcGradient,
                                 anima::TensorCompartment::Pointer &tensorCompartment,
                                 const std::vector <double> &observedData, std::vector <double> &predictedValues);

private:
    DTIEstimationImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    std::vector <double> m_BValuesList;
    std::vector< vnl_vector_fixed<double,3> > m_GradientDirections;

    //! Regular, fast DTI estimation, used to compute starting point for better tensor estimation
    vnl_matrix <double> m_InitialResolutionMatrix;

    double m_B0Threshold;
    typename OutputB0ImageType::Pointer m_EstimatedB0Image;

    static const unsigned int m_NumberOfComponents = 6;
};

} // end of namespace anima

#include "animaDTIEstimationImageFilter.hxx"
