#pragma once

#include <itkImageToImageFilter.h>
#include <vector>
#include <string>

#include <random>

namespace anima
{

template <class TInputImage, class TOutputImage>
class QMRISampleCreationImageFilter
: public itk::ImageToImageFilter <TInputImage, TOutputImage>
{
public:
    /** Standard class typedefs. */
    typedef QMRISampleCreationImageFilter<TInputImage, TOutputImage> Self;
    typedef TInputImage InputImageType;
    typedef TOutputImage OutputImageType;

    typedef itk::ImageToImageFilter<TInputImage,TOutputImage> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    typedef itk::Image <unsigned short, 3> MaskImageType;
    typedef MaskImageType::Pointer MaskImagePointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Type macro that defines a name for this class. */
    itkTypeMacro(QMRISampleCreationImageFilter, ImageToImageFilter)

    /** Smart pointer typedef support.  */
    typedef typename TInputImage::Pointer       InputImagePointer;
    typedef typename TInputImage::ConstPointer  InputImageConstPointer;

    itkSetMacro(QMRILesionRelationshipsFile, std::string)
    void ReadLesionSizesDistributions(std::string sizeDistributionFile);

    void AddQMRIVarianceImage(TInputImage *varImage);

    void SetLesionsProbabilityMap (TInputImage *probaImage)
    {
        m_LesionsProbabilityMap = probaImage;
    }

    itkSetMacro(LesionDiffusionThreshold, double)
    itkSetMacro(LesionMinimalSize, unsigned int)
    itkSetMacro(MinimalDistanceBetweenLesions, double)
    itkSetMacro(NumberOfSeeds, unsigned int)

    itkGetMacro(LesionsOutputMask, MaskImageType *)

protected:
    QMRISampleCreationImageFilter();
    virtual ~QMRISampleCreationImageFilter() {}

    void GenerateData() ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(QMRISampleCreationImageFilter);

    void CheckDataCoherence();
    void InitializeOutputs();
    void ReadQMRILesionRelationships();

    void GenerateQMRIHealthySamples();
    void GenerateAndGrowLesions();
    void UpdateQMRIOnLesions();
    double GetRandomLesionSizeFromDistribution();

    std::vector <InputImagePointer> m_QMRIStdevImages;
    InputImagePointer m_LesionsProbabilityMap;

    std::vector <double> m_XAxisLesionSizesDistribution;
    std::vector <double> m_YAxisLesionSizesDistribution;

    MaskImagePointer m_LesionsOutputMask;

    //! Linear relationship between lesion size and number of iterations in diffusion: y=Ax + B
    static const double m_LesionSizeAFactor, m_LesionSizeBFactor;

    //! Threshold for diffused lesions
    double m_LesionDiffusionThreshold;

    double m_MinimalDistanceBetweenLesions;
    unsigned int m_LesionMinimalSize;

    //! Gaussian relationship between qMRI values inside and outside lesion: I / O ~ N(a,b^2)
    std::string m_QMRILesionRelationshipsFile;
    std::vector <double> m_QMRILesionMeanRelationships;
    vnl_matrix <double> m_QMRILesionCovarianceRelationship;

    unsigned int m_NumberOfSeeds;
    unsigned int m_NumberOfIndividualLesionsKept;

    //! Random generator
    std::mt19937 m_Generator;
};

} // end of namespace anima

#include "animaQMRISampleCreationImageFilter.hxx"
