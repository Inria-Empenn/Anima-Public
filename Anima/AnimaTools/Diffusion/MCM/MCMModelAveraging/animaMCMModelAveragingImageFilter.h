#pragma once

#include <animaNumberedThreadImageToImageFilter.h>
#include <animaMCMImage.h>
#include <itkImage.h>
#include <itkSymmetricEigenAnalysis.h>

#include <animaBaseCompartment.h>
#include <animaMultiCompartmentModel.h>

#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix_fixed.h>

namespace anima
{

template <class PixelScalarType>
class MCMModelAveragingImageFilter :
public anima::NumberedThreadImageToImageFilter < anima::MCMImage<PixelScalarType,3>, anima::MCMImage<PixelScalarType,3> >
{
public:
    /** Standard class typedefs. */
    typedef MCMModelAveragingImageFilter Self;
    typedef anima::MCMImage<PixelScalarType,3> InputImageType;
    typedef anima::MCMImage<PixelScalarType,3> OutputImageType;
    typedef itk::Image<PixelScalarType,3> ScalarImageType;
    typedef itk::Image<unsigned int,3> MoseImageType;
    typedef itk::Image<PixelScalarType,4> Image4DType;
    typedef anima::NumberedThreadImageToImageFilter <InputImageType, OutputImageType> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    typedef anima::MultiCompartmentModel MCModelType;
    typedef typename MCModelType::Pointer MCModelPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(MCMModelAveragingImageFilter, anima::NumberedThreadImageToImageFilter)

    /** Image typedef support */
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;
    typedef typename ScalarImageType::Pointer ScalarImagePointer;
    typedef typename OutputImageType::PixelType OutputPixelType;

    /** Superclass typedefs. */
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    /** Tensor typdefs. */
    typedef vnl_matrix_fixed <double,3,3> MatrixType;
    typedef vnl_vector_fixed <double,3> VectorType;
    typedef itk::SymmetricEigenAnalysis <MatrixType,VectorType,MatrixType> EigenAnalysisType;

    typedef anima::BaseCompartment BaseCompartmentType;
    typedef std::vector <BaseCompartmentType> ModelDataType;
    typedef std::vector <double> WeightDataType;

    itkSetMacro(WeightThreshold, double)
    void SetAICcVolume(unsigned int i, ScalarImageType *vol);
    void SetB0Volume(unsigned int i, ScalarImageType *vol);
    void SetNoiseVolume(unsigned int i, ScalarImageType *vol);

    itkSetMacro(SquaredSimilarity,bool)
    itkSetMacro(SimplifyModels,bool)

    MoseImageType *GetMoseMap () {return m_MoseMap;}
    ScalarImageType *GetOutputB0Volume () {return m_OutputB0Volume;}
    ScalarImageType *GetOutputNoiseVolume () {return m_OutputNoiseVolume;}

protected:
    MCMModelAveragingImageFilter()
    : Superclass()
    {
        m_WeightThreshold = 0.05;
        m_SquaredSimilarity = false;
        m_SimplifyModels = false;
    }

    virtual ~MCMModelAveragingImageFilter() {}

    void GenerateOutputInformation() ITK_OVERRIDE;
    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;
    void AfterThreadedGenerateData() ITK_OVERRIDE;

    void InitializeReferenceOutputModel();
    void IncrementModelPairingVector(std::vector <unsigned int> &modelPairingVector);
    WeightDataType GetAkaikeWeights(const WeightDataType &aicData);

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(MCMModelAveragingImageFilter);

    std::vector <MCModelPointer> m_ReferenceModels;
    MCModelPointer m_ReferenceOutputModel;
    std::vector <ScalarImagePointer> m_AICcVolumes;
    std::vector <ScalarImagePointer> m_B0Volumes;
    std::vector <ScalarImagePointer> m_NoiseVolumes;
    std::vector <unsigned int> m_WorkNonFreeWaterCorrespondences, m_NumberOfNonFreeWaterCompartments;    

    unsigned int m_NumberOfIsotropicCompartments;

    ScalarImagePointer m_OutputB0Volume;
    ScalarImagePointer m_OutputNoiseVolume;
    MoseImageType::Pointer m_MoseMap;

    double m_WeightThreshold;
    bool m_SquaredSimilarity;
    bool m_SimplifyModels;

    static const double m_ZeroThreshold;
};
    
} // end of namespace anima

#include "animaMCMModelAveragingImageFilter.hxx"
