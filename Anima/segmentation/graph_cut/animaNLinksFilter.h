#pragma once

#include "animaGraph.h"
#include <itkCSVArray2DFileReader.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkImageToImageFilter.h>
#include <itkVariableSizeMatrix.h>

namespace anima {
/**
 * @brief Class creating a 3D graph in a graph cut framework
 *
 * The nodes correspond to the voxels of the image to segment.
 * Two additional nodes represent the object (SOURCE) and the background (SINK).
 * N-links represent edges between neighboring nodes/voxels.
 * T-links that bind each classical node to both the SOURCE and the SINK
 * represent the probability for the corresponding voxel to belong respectively
 * to the object and to the background.
 *
 */
template <typename TInput, typename TOutput>
class NLinksFilter : public itk::ImageToImageFilter<TInput, TOutput> {
public:
  /** Standard class typedefs. */
  using Self = NLinksFilter;
  using Superclass = itk::ImageToImageFilter<TInput, TOutput>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(NLinksFilter, ImageToImageFilter);

  /** Image typedef support */
  using InputImagePointer = typename TInput::Pointer;
  using InputImageConstPointer = typename TInput::ConstPointer;
  using InRegionConstIteratorType = itk::ImageRegionConstIterator<TInput>;
  using InputPixelType = typename TInput::PixelType;

  using OutputImagePointer = typename TOutput::Pointer;
  using OutputPixelType = typename TOutput::PixelType;
  using OutputIteratorType = typename itk::ImageRegionIterator<TOutput>;

  using PixelTypeD = double;
  using TSeedProba = itk::Image<PixelTypeD, 3>;
  using ImagePointerProba = TSeedProba::Pointer;
  using ImageIteratorTypeProba = itk::ImageRegionIterator<TSeedProba>;
  using SeedProbaPixelType = TSeedProba::PixelType;

  using PixelTypeInt = int;
  using ImageTypeInt = itk::Image<PixelTypeInt, 3>;
  using pixelIndexInt = ImageTypeInt::IndexType;
  using ImageIteratorTypeInt = itk::ImageRegionIterator<ImageTypeInt>;

  using PixelTypeUC = unsigned char;
  using TMask = itk::Image<PixelTypeUC, 3>;
  using TMaskPointer = TMask::Pointer;
  using MaskRegionIteratorType = itk::ImageRegionIterator<TMask>;
  using MaskRegionConstIteratorType = itk::ImageRegionConstIterator<TMask>;

  using GraphType = Graph<double, double, double>;

  using NumericType = double;
  using doubleVariableSizeMatrixType = itk::VariableSizeMatrix<NumericType>;

  /** The mri images.*/
  void SetInputImage1(const TInput *image);
  void SetInputImage2(const TInput *image);
  void SetInputImage3(const TInput *image);
  void SetInputImage4(const TInput *image);
  void SetInputImage5(const TInput *image);

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

  itkSetMacro(NbModalities, unsigned int);
  itkGetMacro(NbModalities, unsigned int);

  itkSetMacro(Verbose, bool);
  itkGetMacro(Verbose, bool);

protected:
  NLinksFilter() {
    this->SetNthOutput(0, this->MakeOutput(0));
    this->SetNthOutput(1, this->MakeOutput(1));

    m_Sigma = 0.6;
    m_UseSpectralGradient = true;
    m_Verbose = false;
    m_NbModalities = 0;
    m_NbInputs = 4;
    m_NbMaxImage = 10;
    m_IndexImage1 = m_NbMaxImage, m_IndexImage2 = m_NbMaxImage,
    m_IndexImage3 = m_NbMaxImage, m_IndexImage4 = m_NbMaxImage,
    m_IndexImage5 = m_NbMaxImage;

    m_Tol = 0.0001;

    this->SetNumberOfRequiredOutputs(2);
    this->SetNumberOfRequiredInputs(4);
  }

  virtual ~NLinksFilter() {}

  /**  Create the Output */
  itk::DataObject::Pointer MakeOutput(unsigned int idx);

  typename TInput::ConstPointer GetInputImage1();
  typename TInput::ConstPointer GetInputImage2();
  typename TInput::ConstPointer GetInputImage3();
  typename TInput::ConstPointer GetInputImage4();
  typename TInput::ConstPointer GetInputImage5();
  typename TInput::ConstPointer GetInputImage6();

  TMask::ConstPointer GetMask();

  TSeedProba::ConstPointer GetInputSeedProbaSources();
  TSeedProba::ConstPointer GetInputSeedProbaSinks();

  doubleVariableSizeMatrixType GetMatrix(void) { return m_Matrix; }
  std::string GetMatFilename(void) { return m_MatFilename; }

  bool readMatrixFile();

  void CheckSpectralGradient(void);
  void GenerateData() ITK_OVERRIDE;
  void SetGraph();
  bool isInside(unsigned int x, unsigned int y, unsigned int z) const;
  void CreateGraph();
  double computeNLink(int i1, int j1, int k1, int i2, int j2, int k2);

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(NLinksFilter);

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

  typename TOutput::SizeType m_size;

  /** input images: T1 T2 PD FLAIR etc
   */
  std::vector<InputImageConstPointer> m_ListImages;

  bool m_Verbose;

  ImageTypeInt::Pointer pix;

  /** spectral derivatives (e,el,ell)
   */
  TSeedProba::Pointer m_e1, m_e2; // Precomputed spectral grad quantities (keep
                                  // track of 2 images instead of 3...
  unsigned int m_NbModalities;
  unsigned int m_NbInputs, m_NbMaxImage;
  unsigned int m_IndexImage1, m_IndexImage2, m_IndexImage3, m_IndexImage4,
      m_IndexImage5;

  double m_Tol;
};

} // end of namespace anima

#include "animaNLinksFilter.hxx"
