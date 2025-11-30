#pragma once

#include <animaMCMImage.h>
#include <animaNumberedThreadImageToImageFilter.h>
#include <itkImage.h>
#include <itkSymmetricEigenAnalysis.h>

#include <animaBaseCompartment.h>
#include <animaMultiCompartmentModel.h>

#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>

namespace anima {

template <class PixelScalarType>
class MCMModelAveragingImageFilter
    : public anima::NumberedThreadImageToImageFilter<
          anima::MCMImage<PixelScalarType, 3>,
          anima::MCMImage<PixelScalarType, 3>> {
public:
  /** Standard class typedefs. */
  using Self = MCMModelAveragingImageFilter;
  using InputImageType = anima::MCMImage<PixelScalarType, 3>;
  using OutputImageType = anima::MCMImage<PixelScalarType, 3>;
  using ScalarImageType = itk::Image<PixelScalarType, 3>;
  using MoseImageType = itk::Image<unsigned int, 3>;
  using Image4DType = itk::Image<PixelScalarType, 4>;
  using Superclass =
      anima::NumberedThreadImageToImageFilter<InputImageType, OutputImageType>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  using MCModelType = anima::MultiCompartmentModel;
  using MCModelPointer = typename MCModelType::Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(MCMModelAveragingImageFilter,
               anima::NumberedThreadImageToImageFilter);

  /** Image typedef support */
  using InputImagePointer = typename InputImageType::Pointer;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using ScalarImagePointer = typename ScalarImageType::Pointer;
  using OutputPixelType = typename OutputImageType::PixelType;

  /** Superclass typedefs. */
  using InputImageRegionType = typename Superclass::InputImageRegionType;
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  /** Tensor typdefs. */
  using MatrixType = vnl_matrix_fixed<double, 3, 3>;
  using VectorType = vnl_vector_fixed<double, 3>;
  using EigenAnalysisType =
      itk::SymmetricEigenAnalysis<MatrixType, VectorType, MatrixType>;

  using BaseCompartmentType = anima::BaseCompartment;
  using ModelDataType = std::vector<BaseCompartmentType>;
  using WeightDataType = std::vector<double>;

  itkSetMacro(WeightThreshold, double);
  void SetAICcVolume(unsigned int i, ScalarImageType *vol);
  void SetB0Volume(unsigned int i, ScalarImageType *vol);
  void SetNoiseVolume(unsigned int i, ScalarImageType *vol);

  itkSetMacro(SquaredSimilarity, bool);
  itkSetMacro(SimplifyModels, bool);

  MoseImageType *GetMoseMap() { return m_MoseMap; }
  ScalarImageType *GetOutputB0Volume() { return m_OutputB0Volume; }
  ScalarImageType *GetOutputNoiseVolume() { return m_OutputNoiseVolume; }

protected:
  MCMModelAveragingImageFilter() : Superclass() {
    m_WeightThreshold = 0.05;
    m_SquaredSimilarity = false;
    m_SimplifyModels = false;
  }

  virtual ~MCMModelAveragingImageFilter() {}

  void GenerateOutputInformation() ITK_OVERRIDE;
  void BeforeThreadedGenerateData() ITK_OVERRIDE;
  void DynamicThreadedGenerateData(
      const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;
  void AfterThreadedGenerateData() ITK_OVERRIDE;

  void InitializeReferenceOutputModel();
  void
  IncrementModelPairingVector(std::vector<unsigned int> &modelPairingVector);
  WeightDataType GetAkaikeWeights(const WeightDataType &aicData);

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(MCMModelAveragingImageFilter);

  std::vector<MCModelPointer> m_ReferenceModels;
  MCModelPointer m_ReferenceOutputModel;
  std::vector<ScalarImagePointer> m_AICcVolumes;
  std::vector<ScalarImagePointer> m_B0Volumes;
  std::vector<ScalarImagePointer> m_NoiseVolumes;
  std::vector<unsigned int> m_WorkNonFreeWaterCorrespondences,
      m_NumberOfNonFreeWaterCompartments;

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
