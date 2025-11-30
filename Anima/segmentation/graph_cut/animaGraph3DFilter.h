#pragma once

#include "animaNLinksFilter.h"

#include <itkBSplineInterpolateImageFunction.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkIdentityTransform.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkImageToImageFilter.h>
#include <itkMaskImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkVariableSizeMatrix.h>

namespace anima {
/**
 * @brief Class allowing the decimation of the images if necessary (if 3D graph
 * size causes memory problems). This class just launchs NLinksFilter with
 * appropriate image sizes.
 */
template <typename TInput, typename TOutput>
class Graph3DFilter : public itk::ImageToImageFilter<TInput, TOutput> {
public:
  /** Standard class typedefs. */
  using Self = Graph3DFilter;
  using Superclass = itk::ImageToImageFilter<TInput, TOutput>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(Graph3DFilter, ImageToImageFilter);

  /** Image typedef support */
  using InputImagePointer = typename TInput::Pointer;
  using InputImageConstPointer = typename TInput::ConstPointer;
  using InRegionConstIteratorType = itk::ImageRegionConstIterator<TInput>;
  using InRegionIteratorType = itk::ImageRegionIterator<TInput>;
  using InputPixelType = typename TInput::PixelType;

  using PixelTypeD = double;
  using TSeedProba = itk::Image<PixelTypeD, 3>;
  using ImagePointerProba = TSeedProba::Pointer;
  using ImageIteratorTypeProba = itk::ImageRegionIterator<TSeedProba>;
  using SeedProbaPixelType = TSeedProba::PixelType;

  using OutputImagePointer = typename TOutput::Pointer;
  using OutputPixelType = typename TOutput::PixelType;
  using OutputIteratorType = typename itk::ImageRegionIterator<TOutput>;

  using PixelTypeUC = unsigned char;
  using TMask = itk::Image<PixelTypeUC, 3>;
  using TMaskPointer = TMask::Pointer;
  using MaskRegionIteratorType = itk::ImageRegionIterator<TMask>;
  using MaskRegionConstIteratorType = itk::ImageRegionConstIterator<TMask>;

  using GraphType = Graph<double, double, double>;

  using NumericType = double;
  using doubleVariableSizeMatrixType = itk::VariableSizeMatrix<NumericType>;

  using TransformType = itk::IdentityTransform<double, 3>;
  using ResampleImageFilterType = itk::ResampleImageFilter<TInput, TInput>;
  using ResampleImageFilterMaskType = itk::ResampleImageFilter<TMask, TMask>;
  using ResampleImageFilterdoubleType =
      itk::ResampleImageFilter<TSeedProba, TSeedProba>;

  /** The mri images.*/
  void SetInputImage(unsigned int i, const TInput *image);

  /** mask in which the segmentation will be performed
   */
  void SetMask(const TMask *mask);

  /** probabilities of the object and the background, respectively
   */
  void SetInputSeedProbaSources(const TSeedProba *mask);
  void SetInputSeedProbaSinks(const TSeedProba *mask);

  void SetMatFilename(std::string mat) { m_MatFilename = mat; }
  void SetMatrix(doubleVariableSizeMatrixType mat) { m_Matrix = mat; }

  OutputImagePointer GetOutput();
  OutputImagePointer GetOutputBackground();

  TMask::Pointer Upsample(const TMask *input, TMask::SizeType upSize,
                          TMask::DirectionType outputDirection,
                          TMask::PointType outputOrigin);
  TMask::Pointer Dilate(const TMask *input, const TMask *mask);

  /** Superclass typedefs. */
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  void SetTol(const double tol) {
    this->SetCoordinateTolerance(tol);
    this->SetDirectionTolerance(tol);
    m_Tol = tol;
  }

  itkSetMacro(Sigma, double);
  itkGetMacro(Sigma, double);

  itkSetMacro(UseSpectralGradient, bool);
  itkGetMacro(UseSpectralGradient, bool);

  itkSetMacro(Verbose, bool);
  itkGetMacro(Verbose, bool);

protected:
  using NLinksFilterType = NLinksFilter<TInput, TMask>;
  typename NLinksFilterType::Pointer m_NLinksFilter;
  typename NLinksFilterType::Pointer m_NLinksFilterDecim;
  typename ResampleImageFilterType::Pointer m_Resample1;
  typename ResampleImageFilterType::Pointer m_Resample2;
  typename ResampleImageFilterType::Pointer m_Resample3;
  typename ResampleImageFilterType::Pointer m_Resample4;
  typename ResampleImageFilterType::Pointer m_Resample5;
  ResampleImageFilterMaskType::Pointer m_ResampleMask;
  ResampleImageFilterdoubleType::Pointer m_ResampleSources;
  ResampleImageFilterdoubleType::Pointer m_ResampleSinks;

  Graph3DFilter() {
    this->SetNthOutput(0, this->MakeOutput(0));
    this->SetNthOutput(1, this->MakeOutput(1));

    m_Sigma = 0.6;
    m_UseSpectralGradient = true;
    m_Verbose = false;
    m_NbInputs = 4;
    m_NbMaxImages = 10;
    m_IndexImage1 = m_NbMaxImages, m_IndexImage2 = m_NbMaxImages,
    m_IndexImage3 = m_NbMaxImages, m_IndexImage4 = m_NbMaxImages,
    m_IndexImage5 = m_NbMaxImages, m_IndexImage6 = m_NbMaxImages;

    m_NLinksFilter = NLinksFilterType::New();

    m_Tol = 0.0001;

    this->SetNumberOfRequiredOutputs(2);
    this->SetNumberOfRequiredInputs(4);
  }

  virtual ~Graph3DFilter() {}

  /**  Create the Output */
  itk::DataObject::Pointer MakeOutput(unsigned int idx);

  typename TInput::ConstPointer GetInputImage1();
  typename TInput::ConstPointer GetInputImage2();
  typename TInput::ConstPointer GetInputImage3();
  typename TInput::ConstPointer GetInputImage4();
  typename TInput::ConstPointer GetInputImage5();

  TMask::ConstPointer GetMask();

  TSeedProba::ConstPointer GetInputSeedProbaSources();
  TSeedProba::ConstPointer GetInputSeedProbaSinks();

  doubleVariableSizeMatrixType GetMatrix(void) { return m_Matrix; }
  std::string GetMatFilename(void) { return m_MatFilename; }

  void GenerateData() ITK_OVERRIDE;
  bool CheckMemory();
  void ProcessGraphCut();
  void FindDownsampleFactor();
  void InitResampleFilters();
  void ProcessDownsampledGraphCut(int current_count);

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(Graph3DFilter);

  /** set to true to use the spectral gradient instead of a simple gradient
   */
  bool m_UseSpectralGradient;

  /** width of the Gaussian kernel to compute n-links
   */
  double m_Sigma;

  /** the created graph
   */
  GraphType *m_graph;

  /** transformation matrix (from im1,im2,im3 to e,el,ell)
   */
  std::string m_MatFilename;
  doubleVariableSizeMatrixType m_Matrix;

  bool m_Verbose;

  unsigned int m_NbInputs;
  unsigned int m_IndexImage1, m_IndexImage2, m_IndexImage3, m_IndexImage4,
      m_IndexImage5, m_IndexImage6, m_NbMaxImages;

  double m_Tol;

  double m_DownsamplingFactor;
  int m_Count;

  TMask::Pointer m_CurrentMask;
  typename TInput::SizeType m_InputSize;
  typename TInput::DirectionType m_OutputDirection;
  typename TInput::PointType m_OutputOrigin;
};

} // end of namespace anima

#include "animaGraph3DFilter.hxx"
