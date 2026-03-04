#pragma once

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkImageToImageFilter.h>
#include <itkVariableSizeMatrix.h>

enum TLinkMode {
  singleGaussianTLink = 0,
  stremTLink,
};

namespace anima {
/**
 * @brief Class computing the probability maps that are used to create the
 * t-links
 * - single gaussian method computes probability maps using binary masks as
 * input.
 *
 */
template <typename TInput, typename TOutput>
class TLinksFilter : public itk::ImageToImageFilter<TInput, TOutput> {
public:
  /** Standard class typedefs. */
  using Self = TLinksFilter;
  using Superclass = itk::ImageToImageFilter<TInput, TOutput>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(TLinksFilter, ImageToImageFilter);

  /** Image typedef support */
  using InputPixelType = typename TInput::PixelType;
  using InputImagePointer = typename TInput::Pointer;
  using InputImageConstPointer = typename TInput::ConstPointer;
  using InConstIteratorType = itk::ImageRegionConstIterator<TInput>;
  using InIteratorType = itk::ImageRegionIterator<TInput>;

  using OutputPixelType = typename TOutput::PixelType;
  using OutputImagePointer = typename TOutput::Pointer;
  using OutRegionIteratorType = itk::ImageRegionIterator<TOutput>;

  using PixelTypeD = double;
  using TSeedProba = itk::Image<PixelTypeD, 3>;
  using SeedProbaRegionConstIteratorType =
      itk::ImageRegionConstIterator<TSeedProba>;

  using PixelTypeUC = unsigned char;
  using TSeedMask = itk::Image<PixelTypeUC, 3>;
  using TSeedMaskPointer = TSeedMask::Pointer;
  using SeedMaskRegionIteratorType = itk::ImageRegionIterator<TSeedMask>;
  using SeedMaskRegionConstIteratorType =
      itk::ImageRegionConstIterator<TSeedMask>;

  using NumericType = double;
  using DoubleVariableSizeMatrixType = itk::VariableSizeMatrix<NumericType>;
  using MatrixTypeRes = itk::Matrix<double, 1, 1>;

  /** The mri images.*/
  void SetInputImage(unsigned int i, const TInput *image);

  /** probabilities containing the seeds for the object and the background or
   * binary masks containing the seeds for the object and the background
   */
  void SetInputSeedSourcesMask(const TSeedMask *mask);
  void SetInputSeedSinksMask(const TSeedMask *mask);

  void SetInputSeedSourcesProba(const TSeedProba *mask);
  void SetInputSeedSinksProba(const TSeedProba *mask);

  TLinkMode GetTLinkMode() { return m_TLinkMode; }
  void SetTLinkMode(TLinkMode m) { m_TLinkMode = m; }

  TOutput *GetOutputSources();
  TOutput *GetOutputSinks();

  void SetTol(const double tol) {
    this->SetCoordinateTolerance(tol);
    this->SetDirectionTolerance(tol);
  }

  /** Superclass typedefs. */
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  itkSetMacro(Alpha, double);
  itkGetMacro(Alpha, double);

  itkSetMacro(MultiVarSources, double);
  itkGetMacro(MultiVarSources, double);

  itkSetMacro(MultiVarSinks, double);
  itkGetMacro(MultiVarSinks, double);

  itkSetMacro(NbModalities, unsigned int);
  itkGetMacro(NbModalities, unsigned int);

  itkSetMacro(Verbose, bool);
  itkGetMacro(Verbose, bool);

protected:
  TLinksFilter() {
    m_Alpha = 10;
    m_MultiVarSources = 1;
    m_MultiVarSinks = 1;
    m_TLinkMode = singleGaussianTLink;
    m_NbModalities = 0;
    m_NbInputs = 1;
    m_Verbose = false;

    m_NbMaxImage = 11;
    m_IndexSourcesMask = m_NbMaxImage, m_IndexSourcesProba = m_NbMaxImage,
    m_IndexSinksMask = m_NbMaxImage, m_IndexSinksProba = m_NbMaxImage;

    this->SetNumberOfRequiredOutputs(2);
    this->SetNumberOfRequiredInputs(2);

    this->SetNthOutput(0, this->MakeOutput(0));
    this->SetNthOutput(1, this->MakeOutput(1));
  }

  virtual ~TLinksFilter() {}

  TSeedMask::ConstPointer GetInputSeedSourcesMask();
  TSeedMask::ConstPointer GetInputSeedSinksMask();

  TSeedProba::ConstPointer GetInputSeedSourcesProba();
  TSeedProba::ConstPointer GetInputSeedSinksProba();

  void GenerateData() ITK_OVERRIDE;
  void computeSingleGaussian();
  void
  computeSingleGaussianSeeds(TSeedMask::ConstPointer seedMask,
                             OutputImagePointer output, double multiVar,
                             TSeedMask::ConstPointer seedMaskOpp = ITK_NULLPTR);
  void computeStrem();

  /**  Create the Output */
  itk::DataObject::Pointer MakeOutput(unsigned int idx);

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(TLinksFilter);

  /** mixing energies parameters: E_region = m_Alpha * E_intensity
   */
  double m_Alpha;
  double m_MultiVarSources;
  double m_MultiVarSinks;
  TLinkMode m_TLinkMode;
  bool m_Verbose;
  unsigned int m_NbModalities;
  unsigned int m_NbInputs, m_NbMaxImage;
  unsigned int m_IndexSourcesMask, m_IndexSourcesProba, m_IndexSinksMask,
      m_IndexSinksProba;

  std::vector<InConstIteratorType> m_imagesVectorIt;
  std::vector<InIteratorType> m_imagesVectorIt2;
  std::vector<InputImageConstPointer> m_imagesVector;
  std::vector<InputImagePointer> m_imagesVector2;
};

} // end of namespace anima

#include "animaTLinksFilter.hxx"
