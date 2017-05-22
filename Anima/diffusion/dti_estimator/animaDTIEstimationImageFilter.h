#pragma once

#include <animaMaskedImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>

namespace anima
{

template <class InputPixelScalarType, class OutputPixelScalarType>
class DTIEstimationImageFilter :
public anima::MaskedImageToImageFilter < itk::Image<InputPixelScalarType,3>, itk::VectorImage<OutputPixelScalarType,3> >
{
public:
    /** Standard class typedefs. */
    typedef DTIEstimationImageFilter<InputPixelScalarType, OutputPixelScalarType> Self;
    typedef itk::Image<InputPixelScalarType,3> InputImageType;
    typedef itk::VectorImage<OutputPixelScalarType,3> OutputImageType;
    typedef itk::Image<OutputPixelScalarType,3> OutputB0ImageType;
    typedef OutputImageType DTIImageType;
    typedef itk::Image<InputPixelScalarType,4> Image4DType;
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
        vnl_matrix <double> rotationMatrix, workTensor;
        vnl_diag_matrix <double> workEigenValues;
    };

    void SetBValuesList(std::vector <double> bValuesList ) {m_BValuesList = bValuesList;}

    itkSetMacro(B0Threshold, double)
    itkGetMacro(B0Threshold, double)

    itkGetMacro(EstimatedB0Image, OutputB0ImageType *)
    itkGetMacro(EstimatedVarianceImage, OutputB0ImageType *)

    void AddGradientDirection(unsigned int i, vnl_vector_fixed<double,3> &grad);

protected:
    DTIEstimationImageFilter()
        : Superclass()
    {
        m_BValuesList.clear();

        m_B0Threshold = 0;
        m_EstimatedB0Image = NULL;
        m_EstimatedVarianceImage = NULL;
    }

    virtual ~DTIEstimationImageFilter() {}

    void CheckComputationMask() ITK_OVERRIDE;

    void GenerateOutputInformation() ITK_OVERRIDE;
    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId) ITK_OVERRIDE;

    static double OptimizationFunction(const std::vector<double> &x, std::vector<double> &grad, void *func_data);
    double ComputeCostAtPosition(const std::vector<double> &x, const std::vector <double> &observedData,
                                 std::vector <double> &predictedValues, vnl_matrix <double> &rotationMatrix,
                                 vnl_matrix <double> &workTensor, vnl_diag_matrix <double> &workEigenValues);

    double ComputeB0AndVarianceFromTensorVector(const vnl_matrix <double> &tensorValue, const std::vector <double> &dwiSignal, double &outVarianceValue);

private:
    DTIEstimationImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    std::vector <double> m_BValuesList;
    std::vector< vnl_vector_fixed<double,3> > m_GradientDirections;

    double m_B0Threshold;
    typename OutputB0ImageType::Pointer m_EstimatedB0Image, m_EstimatedVarianceImage;

    static const unsigned int m_NumberOfComponents = 6;

    vnl_matrix <double> m_InitialMatrixSolver;
};

} // end of namespace anima

#include "animaDTIEstimationImageFilter.hxx"
