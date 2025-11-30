#pragma once

#include <animaReadWriteFunctions.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkImageToImageFilter.h>
#include <itkLabelContourImageFilter.h>
#include <itkLabelImageToLabelMapFilter.h>
#include <itkLabelMapToLabelImageFilter.h>
#include <itkMultiThreaderBase.h>

namespace anima {
/** @brief Class selecting the connected components touching a given mask
 * border. In MRI, external CSF may contain artifacts due to fluid flow. These
 * effects can cause voxels in the cortex or external CSF to have intensities
 * similar to MS lesions. In order to reduce the number of false positives due
 * to these effects, we remove all candidate lesions that are contiguous to the
 * brain mask border.
 */
template <typename TInput, typename TMask, typename TOutput>
class RemoveTouchingBorderFilter
    : public itk::ImageToImageFilter<TInput, TOutput> {
public:
  /** Standard class typedefs. */
  using Self = RemoveTouchingBorderFilter;
  using Superclass = itk::ImageToImageFilter<TInput, TOutput>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(RemoveTouchingBorderFilter, ImageToImageFilter);

  /** Image typedef support */

  /**  Type of the input image. */
  using InputImageType = TInput;
  using InputImagePointer = typename InputImageType::Pointer;
  using InputIteratorType = typename itk::ImageRegionIterator<InputImageType>;
  using InputConstIteratorType =
      typename itk::ImageRegionConstIterator<InputImageType>;

  /**  Type of the mask image. */
  using MaskImageType = TMask;
  using MaskImagePointer = typename MaskImageType::Pointer;
  using MaskIteratorType = typename itk::ImageRegionIterator<MaskImageType>;
  using MaskConstIteratorType =
      typename itk::ImageRegionConstIterator<MaskImageType>;

  /**  Type of the output image. */
  using OutputImageType = TOutput;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using OutputIteratorType = typename itk::ImageRegionIterator<OutputImageType>;

  using PixelTypeInt = unsigned int;
  using ImageTypeInt = itk::Image<PixelTypeInt, 3>;
  using ImageIteratorTypeInt = itk::ImageRegionIterator<ImageTypeInt>;

  using LabelImageToLabelMapFilterType =
      itk::LabelImageToLabelMapFilter<ImageTypeInt>;
  using TOutputMap = typename LabelImageToLabelMapFilterType::OutputImageType;
  using LabelImageToLabelMapFilterType2 =
      itk::LabelImageToLabelMapFilter<TInput, TOutputMap>;

  using LabelMapToLabelImageFilterType =
      itk::LabelMapToLabelImageFilter<TOutputMap, ImageTypeInt>;

  using ConnectedComponentImageFilterType =
      itk::ConnectedComponentImageFilter<InputImageType, ImageTypeInt>;
  using ConnectedComponentImageFilterType2 =
      itk::ConnectedComponentImageFilter<TMask, ImageTypeInt>;

  using LabelContourImageFilterType =
      itk::LabelContourImageFilter<ImageTypeInt, ImageTypeInt>;

  /** The mri images.*/
  void SetInputImageSeg(const TInput *image);

  void SetOutputNonTouchingBorderFilename(std::string fn) {
    m_OutputNonTouchingBorderFilename = fn;
  }
  void SetOutputTouchingBorderFilename(std::string fn) {
    m_OutputTouchingBorderFilename = fn;
  }

  /** mask in which the segmentation will be performed */
  void SetMask(const TMask *mask);

  TOutput *GetOutputTouchingBorder();
  TOutput *GetOutputNonTouchingBorder();

  void WriteOutputs();

  void SetTol(const double tol) {
    this->SetCoordinateTolerance(tol);
    this->SetDirectionTolerance(tol);
    m_Tol = tol;
  }

  /** Superclass typedefs. */
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  itkSetMacro(NoContour, bool);
  itkGetMacro(NoContour, bool);

  itkSetMacro(LabeledImage, bool);
  itkGetMacro(LabeledImage, bool);

  itkSetMacro(Verbose, bool);
  itkGetMacro(Verbose, bool);

protected:
  RemoveTouchingBorderFilter() {
    this->SetNumberOfRequiredOutputs(2);
    this->SetNumberOfRequiredInputs(2);

    this->SetNthOutput(0, this->MakeOutput(0));
    this->SetNthOutput(1, this->MakeOutput(1));

    m_LabeledImage = false;
    m_NoContour = false;
    m_Verbose = false;
    m_Tol = 0.0001;
    this->SetNumberOfWorkUnits(
        itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads());
  }

  virtual ~RemoveTouchingBorderFilter() {}

  /**  Create the Output */
  itk::DataObject::Pointer MakeOutput(unsigned int idx);

  typename TInput::ConstPointer GetInputImageSeg();
  typename TMask::ConstPointer GetMask();

  /** Does the real work. */
  virtual void GenerateData() ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(RemoveTouchingBorderFilter);

  bool m_LabeledImage;
  bool m_NoContour;
  bool m_Verbose;

  std::string m_OutputNonTouchingBorderFilename;
  std::string m_OutputTouchingBorderFilename;

  std::vector<int> m_labelsToRemove;

  double m_Tol;
};

} // end of namespace anima

#include "animaRemoveTouchingBorderFilter.hxx"
