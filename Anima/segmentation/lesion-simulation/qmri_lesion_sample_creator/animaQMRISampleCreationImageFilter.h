#pragma once

#include <animaMultivariateNormalDistribution.h>

#include <itkImageToImageFilter.h>

#include <random>
#include <string>
#include <vector>

namespace anima {

template <class TInputImage, class TOutputImage>
class QMRISampleCreationImageFilter
    : public itk::ImageToImageFilter<TInputImage, TOutputImage> {
public:
  /** Standard class typedefs. */
  using Self = QMRISampleCreationImageFilter<TInputImage, TOutputImage>;
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;

  using Superclass = itk::ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  using MaskImageType = itk::Image<unsigned short, 3>;
  using MaskImagePointer = MaskImageType::Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Type macro that defines a name for this class. */
  itkTypeMacro(QMRISampleCreationImageFilter, ImageToImageFilter);

  /** Smart pointer typedef support.  */
  using InputImagePointer = typename TInputImage::Pointer;
  using InputImageConstPointer = typename TInputImage::ConstPointer;

  itkSetMacro(QMRILesionRelationshipsFile, std::string);

  void ReadLesionSizesDistributions(std::string sizeDistributionFile);

  void AddQMRIVarianceImage(TInputImage *varImage);

  void SetLesionsProbabilityMap(TInputImage *probaImage) {
    m_LesionsProbabilityMap = probaImage;
  }

  itkSetMacro(LesionDiffusionThreshold, double);
  itkSetMacro(LesionMinimalSize, unsigned int);
  itkSetMacro(MinimalDistanceBetweenLesions, double);
  itkSetMacro(NumberOfSeeds, unsigned int);

  itkGetMacro(LesionsOutputMask, MaskImageType *);

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

  std::vector<InputImagePointer> m_QMRIStdevImages;
  InputImagePointer m_LesionsProbabilityMap;

  std::vector<double> m_XAxisLesionSizesDistribution;
  std::vector<double> m_YAxisLesionSizesDistribution;

  MaskImagePointer m_LesionsOutputMask;

  //! Linear relationship between lesion size and number of iterations in
  //! diffusion: y=Ax + B
  static const double m_LesionSizeAFactor, m_LesionSizeBFactor;

  //! Threshold for diffused lesions
  double m_LesionDiffusionThreshold;

  double m_MinimalDistanceBetweenLesions;
  unsigned int m_LesionMinimalSize;

  //! Gaussian relationship between qMRI values inside and outside lesion: I / O
  //! ~ N(a,b^2)
  std::string m_QMRILesionRelationshipsFile;
  std::vector<double> m_QMRILesionMeanRelationships;
  vnl_matrix<double> m_QMRILesionCovarianceRelationship;
  anima::MultivariateNormalDistribution m_QMRINormalDistribution;

  unsigned int m_NumberOfSeeds;
  unsigned int m_NumberOfIndividualLesionsKept;

  //! Random generator
  std::mt19937 m_Generator;
};

} // end of namespace anima

#include "animaQMRISampleCreationImageFilter.hxx"
