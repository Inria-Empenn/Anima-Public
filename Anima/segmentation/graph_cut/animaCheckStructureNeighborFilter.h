#pragma once

#include <animaReadWriteFunctions.h>

#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryImageToLabelMapFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkImageToImageFilter.h>
#include <itkLabelContourImageFilter.h>
#include <itkLabelMapToLabelImageFilter.h>
#include <itkRelabelComponentImageFilter.h>

namespace anima {
/** @brief Class removing lesions that are not sufficiently in the white matter.
 * Intensity rules may not be enough to discard false positives, therefore we
 * also use localization information. Considering that MS lesions are typically
 * located in WM, we remove the detected ones that do not sufficiently achieve
 * this condition. This filter take two entries: a lesion mask and a white
 * matter map.
 *
 */
template <typename TInput, typename TMask, typename TOutput = TInput>
class CheckStructureNeighborFilter
    : public itk::ImageToImageFilter<TInput, TOutput> {
public:
  /** Standard class typedefs. */
  using Self = CheckStructureNeighborFilter;
  using Superclass = itk::ImageToImageFilter<TInput, TOutput>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(CheckStructureNeighborFilter, ImageToImageFilter);

  /** Image typedef support */

  /**  Type of the input image. */
  using InputImageType = TInput;
  using InputImageConstPointer = typename InputImageType::ConstPointer;
  using InputPixelType = typename InputImageType::PixelType;
  using InputConstIteratorType =
      typename itk::ImageRegionConstIterator<InputImageType>;

  /**  Type of the map image. */
  using MaskImageType = TMask;
  using MaskImageConstPointer = typename MaskImageType::ConstPointer;
  using MaskPixelType = typename MaskImageType::PixelType;
  using MaskConstIteratorType =
      typename itk::ImageRegionConstIterator<MaskImageType>;

  /**  Type of the mask image. */
  using OutputImageType = TOutput;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using OutputPixelType = typename OutputImageType::PixelType;
  using OutputIteratorType = typename itk::ImageRegionIterator<OutputImageType>;

  using PixelTypeInt = int;
  using ImageTypeInt = itk::Image<PixelTypeInt, 3>;
  using ImageIteratorTypeInt = itk::ImageRegionIterator<ImageTypeInt>;

  using ConnectedComponentFilterType =
      typename itk::ConnectedComponentImageFilter<InputImageType, ImageTypeInt>;
  using LabelContourFilterType =
      itk::LabelContourImageFilter<ImageTypeInt, ImageTypeInt>;
  using StructuringElementType =
      itk::BinaryBallStructuringElement<PixelTypeInt, 3>;
  using DilateFilterType =
      itk::GrayscaleDilateImageFilter<ImageTypeInt, ImageTypeInt,
                                      StructuringElementType>;

  /** The mri images.*/
  void SetInputClassification(const TInput *image);
  void SetInputMap(const TMask *image);

  void WriteOutputs();

  void SetTol(const double tol) {
    this->SetCoordinateTolerance(tol);
    this->SetDirectionTolerance(tol);
    m_Tol = tol;
  }

  /** Superclass typedefs. */
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  itkSetMacro(LabelToCheck, int);
  itkGetMacro(LabelToCheck, int);

  itkSetMacro(Verbose, bool);
  itkGetMacro(Verbose, bool);

  itkSetMacro(Ratio, double);
  itkGetMacro(Ratio, double);

  std::string GetOutputFilename() { return m_OutputFilename; }
  void SetOutputFilename(std::string filename) { m_OutputFilename = filename; }

protected:
  CheckStructureNeighborFilter() {
    this->SetNumberOfRequiredOutputs(1);
    this->SetNumberOfRequiredInputs(2);

    m_Verbose = false;
    m_Tol = 0.0001;
    m_Ratio = 0;

    this->SetNumberOfWorkUnits(
        itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads());
  }

  virtual ~CheckStructureNeighborFilter() {}

  typename TInput::ConstPointer GetInputClassification();
  typename TMask::ConstPointer GetInputMap();

  void GenerateData() ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(CheckStructureNeighborFilter);

  std::string m_OutputFilename;
  int m_LabelToCheck;
  bool m_Verbose;
  double m_Tol;
  double m_Ratio;
};

} // end of namespace anima

#include "animaCheckStructureNeighborFilter.hxx"
