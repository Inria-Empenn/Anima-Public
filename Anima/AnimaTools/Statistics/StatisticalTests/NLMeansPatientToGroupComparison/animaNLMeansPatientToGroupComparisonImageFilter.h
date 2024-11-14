#pragma once

#include <iostream>
#include <animaMaskedImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>

#include <vector>

namespace anima
{

template <class PixelScalarType>
class NLMeansPatientToGroupComparisonImageFilter :
        public anima::MaskedImageToImageFilter< itk::VectorImage <PixelScalarType, 3> , itk::Image <PixelScalarType, 3> >
{
public:
    /** Standard class typedefs. */
    typedef NLMeansPatientToGroupComparisonImageFilter<PixelScalarType> Self;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(NLMeansPatientToGroupComparisonImageFilter, MaskedImageToImageFilter)

    /** Image typedef support */
    typedef itk::VectorImage <PixelScalarType, 3> InputImageType;
    typedef itk::Image <PixelScalarType, 3> OutputImageType;

    typedef typename InputImageType::PixelType VectorType;

    typedef vnl_matrix <double> CovarianceType;

    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename InputImageType::IndexType InputImageIndexType;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    /** Superclass typedefs. */
    typedef anima::MaskedImageToImageFilter< InputImageType, OutputImageType > Superclass;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
    typedef typename Superclass::MaskImageType MaskImageType;

    typedef struct
    {
        OutputImageRegionType movingRegion;
        unsigned int numDatabaseImage;
        InputImageIndexType centerIndex;
    } SampleIdentifier;

    void AddDatabaseInput(InputImageType *tmpIm)
    {
        m_DatabaseImages.push_back(tmpIm);
    }

    itkSetMacro(DatabaseCovarianceDistanceAverage, OutputImagePointer)
    itkSetMacro(DatabaseCovarianceDistanceStd, OutputImagePointer)
    itkSetMacro(DatabaseMeanDistanceAverage, OutputImagePointer)
    itkSetMacro(DatabaseMeanDistanceStd, OutputImagePointer)

    itkSetMacro(PatchHalfSize, unsigned int)
    itkSetMacro(SearchNeighborhood, unsigned int)
    itkSetMacro(SearchStepSize, unsigned int)
    itkSetMacro(WeightThreshold, double)
    itkSetMacro(MeanThreshold, double)
    itkSetMacro(VarianceThreshold, double)
    itkSetMacro(BetaParameter, double)

protected:
    NLMeansPatientToGroupComparisonImageFilter()
        : Superclass()
    {
        this->SetNumberOfRequiredOutputs(3);
        this->SetNthOutput(0,this->MakeOutput(0));
        this->SetNthOutput(1,this->MakeOutput(1));
        this->SetNthOutput(2,this->MakeOutput(2));

        m_DatabaseImages.clear();

        m_DatabaseCovarianceDistanceAverage = NULL;
        m_DatabaseCovarianceDistanceStd = NULL;
        m_DatabaseMeanDistanceAverage = NULL;
        m_DatabaseMeanDistanceStd = NULL;

        m_WeightThreshold = 0.0;
        m_MeanThreshold = 0.5;
        m_VarianceThreshold = 6.0;

        m_BetaParameter = 1;
        m_PatchHalfSize = 1;
        m_SearchStepSize = 2;
        m_SearchNeighborhood = 4;
    }

    virtual ~NLMeansPatientToGroupComparisonImageFilter() {}

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(NLMeansPatientToGroupComparisonImageFilter);

    double ComputeWeightedDistanceScore(itk::VariableLengthVector <double> &patientSample, std::vector < double > &databaseWeights,
                                        std::vector < itk::VariableLengthVector <double> > &databaseSamples, double &diffScore);

    std::vector <InputImagePointer> m_DatabaseImages;

    OutputImagePointer m_DatabaseCovarianceDistanceAverage;
    OutputImagePointer m_DatabaseCovarianceDistanceStd;
    OutputImagePointer m_DatabaseMeanDistanceAverage;
    OutputImagePointer m_DatabaseMeanDistanceStd;

    double m_WeightThreshold;
    double m_MeanThreshold, m_VarianceThreshold;
    double m_BetaParameter;

    unsigned int m_PatchHalfSize;
    unsigned int m_SearchStepSize;
    unsigned int m_SearchNeighborhood;
};

} // end namespace anima

#include "animaNLMeansPatientToGroupComparisonImageFilter.hxx"
