#pragma once

#include <itkImageToImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkBinaryImageToLabelMapFilter.h>
#include <itkLabelMapToLabelImageFilter.h>
#include <itkLabelContourImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <animaReadWriteFunctions.h>

namespace anima
{
/** @brief Class selecting the connected components touching a given mask border.
 * In MRI, external CSF may contain artifacts due to fluid flow.
 * These effects can cause voxels in the cortex or external CSF to have intensities similar to MS lesions.
 * In order to reduce the number of false positives due to these effects, we remove all candidate lesions
 * that are contiguous to the brain mask border.
 */
    template<typename TInput, typename TMask, typename TOutput = TInput>
	class RemoveTouchingBorderFilter :
        public itk::ImageToImageFilter< TInput , TOutput >
	{
	  public:
	  /** Standard class typedefs. */
	  typedef RemoveTouchingBorderFilter Self;
	  typedef itk::ImageToImageFilter< TInput , TOutput > Superclass;
	  typedef itk::SmartPointer<Self> Pointer;
	  typedef itk::SmartPointer<const Self>  ConstPointer;

	  /** Method for creation through the object factory. */
	  itkNewMacro(Self);

	  /** Run-time type information (and related methods) */
	  itkTypeMacro(RemoveTouchingBorderFilter, ImageToImageFilter);

	  /** Image typedef support */

	  /**  Type of the input image. */
	  typedef TInput InputImageType;
	  typedef typename InputImageType::Pointer InputImagePointer;
	  typedef typename itk::ImageRegionIterator< InputImageType > InputIteratorType;
	  typedef typename itk::ImageRegionConstIterator< InputImageType > InputConstIteratorType;

	  /**  Type of the mask image. */
	  typedef TMask MaskImageType;
	  typedef typename MaskImageType::Pointer MaskImagePointer;
	  typedef typename itk::ImageRegionIterator< MaskImageType > MaskIteratorType;
	  typedef typename itk::ImageRegionConstIterator< MaskImageType > MaskConstIteratorType;

	  /**  Type of the output image. */
	  typedef TOutput OutputImageType;
	  typedef typename OutputImageType::Pointer OutputImagePointer;
	  typedef typename itk::ImageRegionIterator< OutputImageType > OutputIteratorType;

	  typedef int 	PixelTypeInt;
	  typedef itk::Image <PixelTypeInt,3> ImageTypeInt;
	  typedef itk::ImageRegionIterator< ImageTypeInt > ImageIteratorTypeInt;

	  typedef itk::BinaryImageToLabelMapFilter<TInput> BinaryImageToLabelMapFilterType;
	  typedef typename BinaryImageToLabelMapFilterType::OutputImageType TOutputMap;
	  typedef itk::LabelMapToLabelImageFilter<TOutputMap, ImageTypeInt> LabelMapToLabelImageFilterType;

	  typedef itk::LabelMapToLabelImageFilter<TOutputMap, OutputImageType> LabelMapToLabelTOutImageFilterType;

	  /** The mri images.*/
	  void SetInputImageSeg(const TInput* image);

	  void SetOutputNonTouchingBorderFilename(std::string fn){m_OutputNonTouchingBorderFilename = fn;}
	  void SetOutputTouchingBorderFilename(std::string fn){m_OutputTouchingBorderFilename = fn;}

	  /** mask in which the segmentation will be performed
	   */
	  void SetMask(const TMask* mask);

	  TOutput* GetOutputTouchingBorder();
	  TOutput* GetOutputNonTouchingBorder();

	  void WriteOutputs();

	  void SetTol(const double tol)
	  {
	      this->SetCoordinateTolerance(tol);
	      this->SetDirectionTolerance(tol);
	      m_Tol = tol;
	  }

	  /** Superclass typedefs. */
	  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

	  itkSetMacro(Verbose, bool);
	  itkGetMacro(Verbose, bool);


	  protected:

	  RemoveTouchingBorderFilter()
	  {

	      this->SetNumberOfRequiredOutputs(2);
	      this->SetNumberOfRequiredInputs(2);

	      this->SetNthOutput( 0, this->MakeOutput(0) );
	      this->SetNthOutput( 1, this->MakeOutput(1) );

	      m_Verbose=false;
	      m_Tol = 0.0001;
	      this->SetNumberOfThreads(itk::MultiThreader::GetGlobalDefaultNumberOfThreads());

	  }

	  virtual ~RemoveTouchingBorderFilter()
	  {
	  }

	  /**  Create the Output */
	  itk::DataObject::Pointer MakeOutput(unsigned int idx);

	  typename TInput::ConstPointer GetInputImageSeg();
	  typename TMask::ConstPointer GetMask();

	  /** Does the real work. */
	  virtual void GenerateData() ITK_OVERRIDE;

	  private:
	  RemoveTouchingBorderFilter(const Self&); //purposely not implemented
	  void operator=(const Self&); //purposely not implemented

	  bool m_Verbose;

	  std::string m_OutputNonTouchingBorderFilename;
	  std::string m_OutputTouchingBorderFilename;

	  std::vector<int> m_labelsToRemove;

	  double m_Tol;

	};
} // end of namespace anima

#include "animaRemoveTouchingBorderFilter.hxx"
