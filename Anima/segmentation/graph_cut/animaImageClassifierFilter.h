#pragma once

#include <animaReadWriteFunctions.h>
#include <itkGaussianMembershipFunction.h>
#include <itkImageToImageFilter.h>

namespace anima {
/**
 * @brief Classify each voxels into one of the given GMM classes.
 *
 *
 */
template <typename TInput, typename TMask, typename TOutput = TInput>
class ImageClassifierFilter : public itk::ImageToImageFilter<TInput, TOutput> {
public:
  /** Standard class typedefs. */
  using Self = ImageClassifierFilter;
  using Superclass = itk::ImageToImageFilter<TInput, TOutput>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(ImageClassifierFilter, ImageToImageFilter);

  /** Image typedef support */
  using InputImageConstPointer = typename TInput::ConstPointer;
  using InputConstIteratorType = typename itk::ImageRegionConstIterator<TInput>;

  /** Mask typedef support */
  using MaskConstIteratorType = typename itk::ImageRegionConstIterator<TMask>;

  using OutputImagePointer = typename TOutput::Pointer;
  using OutputImageConstPointer = typename TOutput::ConstPointer;
  using OutputIteratorType = typename itk::ImageRegionIterator<TOutput>;

  using NumericType = double;
  using DoubleVariableSizeMatrixType = itk::VariableSizeMatrix<NumericType>;

  using MeasurementVectorType = itk::VariableLengthVector<double>;
  using GaussianFunctionType =
      itk::Statistics::GaussianMembershipFunction<MeasurementVectorType>;

  void SetTol(const double tol) {
    this->SetCoordinateTolerance(tol);
    this->SetDirectionTolerance(tol);
  }

  /** The mri images.*/
  void SetInputImage1(const TInput *image);
  void SetInputImage2(const TInput *image);
  void SetInputImage3(const TInput *image);
  void SetInputImage4(const TInput *image);
  void SetInputImage5(const TInput *image);

  /** mask in which the segmentation will be performed
   */
  void SetMask(const TMask *mask);

  /** Superclass typedefs. */
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  itkSetMacro(Verbose, bool);
  itkGetMacro(Verbose, bool);

  void WriteOutputs();
  void SetGaussianModel(std::vector<GaussianFunctionType::Pointer> &model) {
    m_GaussianModel = model;
  }
  void SetAlphas(std::vector<double> &model) { m_Alphas = model; }
  void SetOutputFilename(std::string filename) { m_OutputFilename = filename; }
  std::vector<GaussianFunctionType::Pointer> GetGaussiabModel() {
    return m_GaussianModel;
  }
  std::vector<double> GetAlphas() { return m_Alphas; }

protected:
  ImageClassifierFilter() {
    this->SetNumberOfRequiredOutputs(1);
    this->SetNumberOfRequiredInputs(2);

    m_NbInputs = 2;
    m_NbMaxImages = 7;
    m_IndexImage1 = m_NbMaxImages, m_IndexImage2 = m_NbMaxImages,
    m_IndexImage3 = m_NbMaxImages, m_IndexImage4 = m_NbMaxImages,
    m_IndexImage5 = m_NbMaxImages, m_IndexImage6 = m_NbMaxImages;
    m_Verbose = false;
    this->SetNumberOfWorkUnits(
        itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads());
  }

  virtual ~ImageClassifierFilter() {}

  typename TInput::ConstPointer GetInputImage1();
  typename TInput::ConstPointer GetInputImage2();
  typename TInput::ConstPointer GetInputImage3();
  typename TInput::ConstPointer GetInputImage4();
  typename TInput::ConstPointer GetInputImage5();

  typename TMask::ConstPointer GetMask();

  void DynamicThreadedGenerateData(
      const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;
  double probability(DoubleVariableSizeMatrixType &point,
                     GaussianFunctionType::Pointer model);

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(ImageClassifierFilter);

  std::vector<double> m_Alphas;
  std::vector<GaussianFunctionType::Pointer> m_GaussianModel;

  std::vector<typename TInput::ConstPointer> m_ImagesVector;
  std::string m_OutputFilename;

  unsigned int m_NbInputs, m_NbMaxImages;
  unsigned int m_IndexImage1, m_IndexImage2, m_IndexImage3, m_IndexImage4,
      m_IndexImage5, m_IndexImage6;

  bool m_Verbose;
};

} // end of namespace anima

#include "animaImageClassifierFilter.hxx"
